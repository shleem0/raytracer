import bpy
import json
import mathutils
import os

def get_camera_data(camera):
    """Get camera data"""
    location = camera.location
    camera_matrix = camera.matrix_world
    forward = camera_matrix.to_3x3() @ mathutils.Vector((0, 0, -1))
    gaze_vector = {
        "x": forward.x,
        "y": forward.y,
        "z": forward.z
    }
    up_vector_world = camera_matrix.to_3x3() @ mathutils.Vector((0, 1, 0))
    up_vector = {
        "x": up_vector_world.x,
        "y": up_vector_world.y,
        "z": up_vector_world.z
    }
    
    focal_length = camera.data.lens
    
    sensor = {
        "width": camera.data.sensor_width,
        "height": camera.data.sensor_height
    }
    
    scene = bpy.context.scene
    
    film_resolution = {
        "width": scene.render.resolution_x,
        "height": scene.render.resolution_y
    }
    
    return {
        "location": {"x": location.x, "y": location.y, "z": location.z},
        "gaze_vector": gaze_vector,
        "up_vector": up_vector,
        "focal_length": focal_length,
        "sensor": sensor,
        "film_resolution": film_resolution
    }



def get_point_light_data(point_light):
    """Get point light data"""
    location = point_light.location
    radiant_intensity = point_light.data.energy
    return {
        "location": {"x": location.x, "y": location.y, "z": location.z},
        "radiant_intensity": radiant_intensity
    }

def get_sphere_data(sphere):
    """Get sphere data"""
    location = sphere.location
    radius = sphere.scale.x  # spheres are often scaled equally in all directions
    material = get_material_data(sphere)
    return {
        "location": {"x": location.x, "y": location.y, "z": location.z},
        "radius": radius,
        "material": material
    }

def get_cube_data(cube):
    """Get cube data"""
    translation = cube.location
    rotation = cube.rotation_euler
    scale = cube.scale.x  # cubes are often scaled equally in all directions
    material = get_material_data(cube)
    return {
        "translation": {"x": translation.x, "y": translation.y, "z": translation.z},
        "rotation": {"x": rotation.x, "y": rotation.y, "z": rotation.z},
        "scale": scale,
        "material": material
    }

def get_plane_data(plane):
    """Get plane data"""
    # Get the corners of the plane by applying the matrix_world to the vertices
    corners = []
    material = get_material_data(plane)
    for vertex in plane.data.vertices:
        corner = plane.matrix_world @ vertex.co
        corners.append({"x": corner.x, "y": corner.y, "z": corner.z})
    return {
        "corners": corners,
        "material": material
    }


def get_material_data(obj):
    """Extract material info from a mesh object"""
    if not obj.data.materials or obj.data.materials[0] is None:
        # Default material if none assigned
        return {
            "diffuse": {"r": 1.0, "g": 1.0, "b": 1.0},
            "specular": {"r": 1.0, "g": 1.0, "b": 1.0},
            "shininess": 32.0,
            "transparency": 0.0,
            "ior": 1.0,
            "texture": None
        }

    mat = obj.data.materials[0]  # Take the first material for simplicity

    # Default values
    diffuse = {"r": 1.0, "g": 1.0, "b": 1.0}
    specular = {"r": 1.0, "g": 1.0, "b": 1.0}
    shininess = 32.0
    transparency = 0.0
    ior = 1.0
    texture_path = None

    if mat.node_tree:
        for node in mat.node_tree.nodes:
            if node.type == 'BSDF_PRINCIPLED':
                # Diffuse color

                if 'Base Color' in node.inputs:
                    colour_node = node.inputs['Base Color']
                    base_color = colour_node.default_value
                    diffuse = {"r": base_color[0], "g": base_color[1], "b": base_color[2]}
                    texture_path = get_material_texture(mat)

                
                # Specular intensity (scalar)
                if 'Specular' in node.inputs:
                    s = node.inputs['Specular'].default_value
                    specular = {"r": s, "g": s, "b": s}
                
                # Roughness → shininess
                if 'Roughness' in node.inputs:
                    roughness = node.inputs['Roughness'].default_value
                    shininess = max(1.0, (1.0 - roughness) * 128.0)
                    
                #Transparency
                if 'Transmission' in node.inputs:
                    transmission = node.inputs['Transmission'].default_value
                    transparency = max(0.0, min(1.0, transmission))

                # --- Refractive Index (IOR) ---
                if 'IOR' in node.inputs:
                    ior = node.inputs['IOR'].default_value

                break  # Only process first Principled BSDF node

    return {
        "diffuse": diffuse,
        "specular": specular,
        "shininess": shininess,
        "transparency": transparency,
        "ior": ior,
        "texture": texture_path if texture_path else None
    }


def get_material_texture(mat):
    """Return the filepath of the first Image Texture node in the material, if any."""
    if not mat.node_tree:
        return None

    for node in mat.node_tree.nodes:
        if node.type == 'TEX_IMAGE' and node.image:
            # Return absolute path
            return os.path.basename(node.image.filepath_raw) or node.image.name
    return None


def export_to_json():
    """Export the scene data to JSON"""
    cameras = []
    point_lights = []
    spheres = []
    cubes = []
    planes = []

    for object in bpy.context.scene.objects:
        
        print(object.data.name)
        
        if object.type == 'CAMERA':
            cameras.append(get_camera_data(object))
        elif object.type == 'LIGHT' and object.data.type == 'POINT':
            point_lights.append(get_point_light_data(object))
        elif (object.type == 'MESH' and object.data.name.startswith('Sphere')) or object.type == "META":
            spheres.append(get_sphere_data(object))
        elif object.type == 'MESH' and object.data.name.startswith('Cube'):
            cubes.append(get_cube_data(object))
        elif object.type == 'MESH' and object.data.name.startswith('Plane'):
            planes.append(get_plane_data(object))

    data = {
        "properties": {
            "cameras": cameras,
            "point_lights": point_lights,
            "spheres": spheres,
            "cubes": cubes,
            "planes": planes
        }
    }
    
    # Get the path to the .blend file
    blend_file_path = bpy.data.filepath
    # Get the directory of the .blend file
    blend_file_dir = os.path.dirname(blend_file_path)

    # Export to JSON
    with open(os.path.join(blend_file_dir, '..\\ASCII\scene.json'), 'w') as f:
        json.dump(data, f, indent=4)

export_to_json()
