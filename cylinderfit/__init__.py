# ##### BEGIN GPL LICENSE BLOCK #####
#
#  CylinderFit, (c) 2024 Michel Anders (varkenvarken)
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "CylinderFit",
    "author": "Michel Anders (varkenvarken)",
    "version": (0, 0, 20240126131613),
    "blender": (4, 0, 0),
    "location": "Edit mode 3d-view, Add-->CylinderFit",
    "description": "Add a cylinder to the mesh that best fits a collection of selected vertices",
    "warning": "",
    "wiki_url": "",
    "category": "Mesh",
}

import numpy as np
import mathutils

import bpy
import bpy.ops
from bpy.ops import mesh


# force reload of sub modules if already/still present
from importlib import reload

for mod in ("fitting",):
    if mod in locals():
        reload(locals()[mod])

from .fitting import fit, fitRod


def cylinderfit(points):
    result = fit(points)
    return (
        mathutils.Vector(result.centroid),
        mathutils.Vector(result.direction),
        result.radius,
        result.fit,
    )


def rodfit(points, quantile):
    result = fitRod(points, quantile)
    return (
        mathutils.Vector(result.centroid),
        mathutils.Vector(result.direction),
        result.radius,
        result.fit,
    )


class CylinderFit(bpy.types.Operator):
    bl_idname = "mesh.cylinderfit"
    bl_label = "CylinderFit"
    bl_options = {"REGISTER", "UNDO"}

    nverts: bpy.props.IntProperty(
        name="Vertices",
        description="Number of cylinder vertices",
        default=32,
        min=3,
        soft_max=128,
    )

    rod: bpy.props.BoolProperty(
        name="Rod",
        description="Assume selected vertices for uniform density rod",
        default=False,
    )

    quantile: bpy.props.FloatProperty(
        name="Radius weight",
        description="Higher is further out",
        min=0.001,
        max=1.0,
        soft_max=0.99,
        default=0.99,
    )

    def draw(self, context):
        layout = self.layout
        layout.prop(self, "nverts")
        row = layout.row()
        row.prop(self,"rod")
        if self.rod:
            row.prop(self,"quantile")
            
    @classmethod
    def poll(self, context):
        return context.mode == "EDIT_MESH" and context.active_object.type == "MESH"

    def execute(self, context):
        bpy.ops.object.editmode_toggle()
        world_matrix = context.active_object.matrix_world
        me = context.active_object.data
        count = len(me.vertices)
        if count > 0:  # degenerate mesh, but better safe than sorry
            # get the vertex coordinates. We retrieve them flattened
            # so we will have to reshape them ourselves
            shape = (count, 3)
            verts = np.empty(count * 3, dtype=np.float32)
            me.vertices.foreach_get("co", verts)
            verts.shape = shape

            # also get the selected status. No need to reshape because it is 1-D
            selected = np.empty(count, dtype=bool)
            me.vertices.foreach_get("select", selected)

            # we don't want to force the user to apply transformations first, so
            # we convert all vertex coordinates to world space, by first adding a
            # column with all ones, so it will match our 4x4 world matrix and
            # translations will work, and multiply the matrix with each vector
            # and finally removing the column of ones again.
            verts = np.c_[verts, np.ones(len(verts))]
            verts = np.dot(world_matrix, verts.T).T[:, :3]

            if np.count_nonzero(selected) >= 6:
                if self.rod:
                    centroid, direction, radius, error = rodfit(verts[selected], self.quantile)
                else:
                    centroid, direction, radius, error = cylinderfit(verts[selected])
                Z = mathutils.Vector((0, 0, 1))
                quaternion = Z.rotation_difference(mathutils.Vector(direction))
                euler = quaternion.to_euler()
                # bpy.ops.object.editmode_toggle()
                print(f"{centroid=}, {direction=}, {radius=}, {error=}")
                z = np.array(Z)
                z.shape = 3, 1
                projections = np.dot(verts[selected] - centroid, direction).flatten()
                # for p, s in zip(projections, verts[selected]):
                #     print(p, s)
                length = np.max(projections) - np.min(projections)
                print(length, np.max(projections), np.min(projections))
                mesh.primitive_cylinder_add(
                    enter_editmode=True,
                    vertices=self.nverts,
                    depth=length,
                    radius=radius,
                    location=centroid,
                    rotation=euler,  # (quaternion @ obj_rotation).to_euler(),
                )
                # context.active_object.rotation_quaternion = obj_rotation
                # context.view_layer.update()
            else:
                self.report(
                    {"WARNING"},
                    "Need at least 6 selected vertices to fit a cylinder through",
                )
                bpy.ops.object.editmode_toggle()
        return {"FINISHED"}


def menu_func(self, context):
    self.layout.operator(
        CylinderFit.bl_idname, text="Fit cylinder to selected", icon="PLUGIN"
    )


def register():
    bpy.utils.register_class(CylinderFit)
    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)


def unregister():
    bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)
    bpy.utils.unregister_class(CylinderFit)


if __name__ == "__main__":
    register()
