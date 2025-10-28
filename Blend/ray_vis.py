import bpy
import mathutils

def visualize_ray(origin, direction):
    # Normalize the direction vector
    direction = direction.normalized()

    # Create a new mesh and object
    mesh = bpy.data.meshes.new("Ray")
    obj = bpy.data.objects.new("Ray", mesh)

    # Define the mesh vertices
    vertices = [(origin.x, origin.y, origin.z), 
                (origin.x + direction.x * 10, origin.y + direction.y * 10, origin.z + direction.z * 10)]

    # Define the mesh edges
    edges = [(0, 1)]

    # Create the mesh from the vertices and edges
    mesh.from_pydata(vertices, edges, [])
    mesh.update()

    # Add the object to the scene
    bpy.context.collection.objects.link(obj)

# Example usage:
origin = mathutils.Vector((0, 0, 0))  # Ray origin
direction = mathutils.Vector((1, 0, 0))  # Ray direction
