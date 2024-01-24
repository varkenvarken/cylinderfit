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
    "version": (0, 0, 20240119100755),
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

from .fitting import fit


def cylinderfit(points):
    direction, centroid, radius, error = fit(points)
    return mathutils.Vector(centroid), mathutils.Vector(direction), radius, error


class CylinderFit(bpy.types.Operator):
    bl_idname = "mesh.cylinderfit"
    bl_label = "CylinderFit"
    bl_options = {"REGISTER", "UNDO"}

    size: bpy.props.FloatProperty(
        name="Length",
        description="Length of the line segment",
        default=1,
        min=0,
        soft_max=10,
    )

    @classmethod
    def poll(self, context):
        return context.mode == "EDIT_MESH" and context.active_object.type == "MESH"

    def execute(self, context):
        # NOTE scale and rotation must be applied!
        bpy.ops.object.editmode_toggle()
        me = context.active_object.data
        count = len(me.vertices)
        if count > 0:  # degenerate mesh, but better safe than sorry
            shape = (count, 3)
            verts = np.empty(count * 3, dtype=np.float32)
            selected = np.empty(count, dtype=np.bool)
            me.vertices.foreach_get("co", verts)
            me.vertices.foreach_get("select", selected)
            verts.shape = shape
            if np.count_nonzero(selected) >= 6:
                centroid, direction, radius, error = cylinderfit(verts[selected])
                Z = mathutils.Vector((0, 0, 1))
                quaternion = Z.rotation_difference(mathutils.Vector(direction))
                euler = quaternion.to_euler()
                nverts = 32
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
                    vertices=nverts,
                    depth=length,
                    location=centroid,
                    radius=radius,
                    rotation=euler,
                )
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
