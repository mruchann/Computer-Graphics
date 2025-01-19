#include <iostream>
#include <cmath>
#include <cfloat>
#include <thread>
#include "parser.h"
#include "ppm.h"

using namespace parser;
using namespace std;

Vec3f getPixel(int i, int j, const Camera& camera, const Vec3f& topLeft, const Vec3f& u, const float& pixelWidth, const float& pixelHeight) {
    float s_u = (0.5 + j) * pixelWidth;
    float s_v = (0.5 + i) * pixelHeight;
    return topLeft + u * s_u - camera.up * s_v;
}

class Ray {
public:
    Vec3f origin, direction;

    Ray(): origin(), direction() {}

    Ray(Vec3f origin, Vec3f direction) {
        this->origin = origin;
        this->direction = direction;
    }
};

float intersect(const Ray& ray, const Sphere& sphere, const Scene& scene) {
    const Vec3f& center = scene.vertex_data[sphere.center_vertex_id - 1];
    float A = ray.direction.square();
    float B = 2 * (ray.direction * (ray.origin - center));
    float C = (ray.origin - center).square() - sphere.radius * sphere.radius;
    float delta = B * B - 4 * A * C;

    if (delta < 0) {
        return FLT_MAX;
    }
    else if (delta == 0) {
        return -B / (2 * A);
    }

    return (- B - sqrt(delta)) / (2 * A);
}

// TODO: note that { col1, col2, col3 }
inline float determinant(vector<Vec3f> mat) {
    return mat[0].x * (mat[1].y * mat[2].z - mat[2].y * mat[1].z)
           - mat[1].x * (mat[0].y * mat[2].z - mat[2].y * mat[0].z)
           + mat[2].x * (mat[0].y * mat[1].z - mat[1].y * mat[0].z);
}

float intersect(const Ray& ray, const Triangle& triangle, const Scene& scene) {
    if (ray.direction * triangle.normal >= 0) {
        return FLT_MAX;
    }
    Vec3f a = scene.vertex_data[triangle.indices.v0_id - 1];
    Vec3f b = scene.vertex_data[triangle.indices.v1_id - 1];
    Vec3f c = scene.vertex_data[triangle.indices.v2_id - 1];
    vector<Vec3f> A = { (a - b), (a - c), ray.direction };
    Vec3f B(a - ray.origin);

    float detA = determinant(A);

    if (detA == 0) {
        return FLT_MAX;
    }

    float beta = determinant({B, A[1], A[2]}) / detA;
    float gamma = determinant({A[0], B, A[2]}) / detA;

    if (beta + gamma < -1e-6 || beta < -1e-6 || gamma < -1e-6 || beta + gamma > 1 + 1e-6 || beta > 1 + 1e-6 || gamma > 1 + 1e-6) {
        return FLT_MAX;
    }

    float t = determinant({A[0], A[1], B}) / detA;

    return t;
}

Vec3f eval(const Ray& ray, float t) {
    return ray.origin + ray.direction * t;
}

Vec3f getColor(const Ray& ray, const Scene& scene, int recursion_depth) {
    float t_min = FLT_MAX;
    int objectType = -1; // 0 -> sphere, 1 -> triangle
    Sphere hitSphere;
    Triangle hitTriangle;

    for (const Sphere& sphere: scene.spheres) {
        float t = intersect(ray, sphere, scene); // returns FLT_MAX if there's no intersection
        if (t < t_min && t >= -1e-6) {
            t_min = t;
            objectType = 0;
            hitSphere = sphere;
        }
    }

    for (const Triangle& triangle: scene.triangles) {
        float t = intersect(ray, triangle, scene); // returns FLT_MAX if there's no intersection
        if (t < t_min && t >= -1e-6) {
            t_min = t;
            objectType = 1;
            hitTriangle = triangle;
        }
    }
    Vec3f color;
    if (objectType == -1) {
        if (recursion_depth == scene.max_recursion_depth) {
            color.x = scene.background_color.x;
            color.y = scene.background_color.y;
            color.z = scene.background_color.z;
        }
        return color;
    }
    if (objectType == 0) { // sphere
        const Material& material = scene.materials[hitSphere.material_id - 1];
        Vec3f P = eval(ray, t_min);
        const Vec3f& center = scene.vertex_data[hitSphere.center_vertex_id - 1];
        Vec3f N = (P - center).normalize();

        // Ambient
        color.x = material.ambient.x * scene.ambient_light.x;
        color.y = material.ambient.y * scene.ambient_light.y;
        color.z = material.ambient.z * scene.ambient_light.z;

        // Reflection
        if (recursion_depth > 0 && material.is_mirror) {
            Vec3f new_direction = (ray.direction + ((N * 2) * max(0.0f, N * -ray.direction))).normalize();
            Ray reflectedRay = Ray(P + N * scene.shadow_ray_epsilon, new_direction);
            Vec3f reflection_color = getColor(reflectedRay, scene, recursion_depth - 1);
            color.x += reflection_color.x * material.mirror.x;
            color.y += reflection_color.y * material.mirror.y;
            color.z += reflection_color.z * material.mirror.z;
        }

        for (const PointLight& light : scene.point_lights) {
            Vec3f S = (light - P);
            bool shadow = false;
            // Shadow
            Ray shadowRay = Ray(P + N * scene.shadow_ray_epsilon, S);
            for (const Sphere& sphere: scene.spheres) {
                float t_shadow = intersect(shadowRay, sphere, scene);
                if (-1e-6 < t_shadow && t_shadow < 1 + 1e-6) {
                    shadow = true;
                    break;
                }
            }
            if (shadow) {
                continue;
            }
            for (const Triangle& triangle: scene.triangles) {
                float t_shadow = intersect(shadowRay, triangle, scene);
                if (-1e-6 < t_shadow && t_shadow < 1 + 1e-6) {
                    shadow = true;
                    break;
                }
            }
            if (shadow) {
                continue;
            }

            Vec3f L = (light - P).normalize();
            Vec3f H = (L - ray.direction).normalize();
            float squared = S * S;
            Vec3f irradiance = light.intensity / squared;
            const float NL = N * L;
            const float NH = N * H;
            // Diffuse
            if (NL > 0) {
                color.x += material.diffuse.x * NL * irradiance.x;
                color.y += material.diffuse.y * NL * irradiance.y;
                color.z += material.diffuse.z * NL * irradiance.z;
            }

            // Specular
            if (NH > 0) {
                const float phong = pow(NH, material.phong_exponent);
                color.x += material.specular.x * phong * irradiance.x;
                color.y += material.specular.y * phong * irradiance.y;
                color.z += material.specular.z * phong * irradiance.z;
            }
        }
        return color;
    }

    // triangle
    const Material& material = scene.materials[hitTriangle.material_id - 1];
    Vec3f P = eval(ray, t_min);
    const Vec3f& N = hitTriangle.normal;

    // Ambient
    color.x = material.ambient.x * scene.ambient_light.x;
    color.y = material.ambient.y * scene.ambient_light.y;
    color.z = material.ambient.z * scene.ambient_light.z;

    // Reflection
    if (recursion_depth > 0 && material.is_mirror) {
        Vec3f new_direction = (ray.direction + ((N * 2) * max(0.0f, N * -ray.direction))).normalize();
        Ray reflectedRay = Ray(P + N * scene.shadow_ray_epsilon, new_direction);
        Vec3f reflection_color = getColor(reflectedRay, scene, recursion_depth - 1);
        color.x += reflection_color.x * material.mirror.x;
        color.y += reflection_color.y * material.mirror.y;
        color.z += reflection_color.z * material.mirror.z;
    }

    for (const PointLight& light: scene.point_lights) {
        Vec3f S = (light - P);
        bool shadow = false;
        // Shadow
        Ray shadowRay = Ray(P + N * scene.shadow_ray_epsilon, S);
        for (const Sphere& sphere: scene.spheres) {
            float t_shadow = intersect(shadowRay, sphere, scene);
            if (-1e-6 < t_shadow && t_shadow < 1 + 1e-6) {
                shadow = true;
                break;
            }
        }
        if (shadow) {
            continue;
        }
        for (const Triangle& triangle: scene.triangles) {
            float t_shadow = intersect(shadowRay, triangle, scene);
            if (-1e-6 < t_shadow && t_shadow < 1 + 1e-6) {
                shadow = true;
                break;
            }
        }
        if (shadow) {
            continue;
        }

        Vec3f L = (light - P).normalize();
        Vec3f H = (L - ray.direction).normalize();
        float squared = S * S;
        Vec3f irradiance = light.intensity / squared;
        const float NL = N * L;
        const float NH = N * H;
        // Diffuse
        if (NL > 0) {
            color.x += material.diffuse.x * NL * irradiance.x;
            color.y += material.diffuse.y * NL * irradiance.y;
            color.z += material.diffuse.z * NL * irradiance.z;
        }

        // Specular
        if (NH > 0) {
            const float phong = pow(NH, material.phong_exponent);
            color.x += material.specular.x * phong * irradiance.x;
            color.y += material.specular.y * phong * irradiance.y;
            color.z += material.specular.z * phong * irradiance.z;
        }
    }
    return color;
}

void work(int l_height, int r_height, int width, const Camera& camera, const Vec3f& topLeft,
          const Vec3f& u, float pixelWidth, float pixelHeight, const Scene& scene, unsigned char* image) {
    for (int y = l_height; y < r_height; ++y) {
        for (int x = 0; x < width; ++x) {
            Vec3f pixel = getPixel(y, x, camera, topLeft, u, pixelWidth, pixelHeight); // s
            Vec3f direction = (pixel - camera.position).normalize(); // d
            Ray ray = Ray(camera.position, direction);
            Vec3f rgb = getColor(ray, scene, scene.max_recursion_depth);

            int ind = (x + y * width) * 3;

            image[ind] = min(255, int(round(rgb.x)));
            image[ind+1] = min(255, int(round(rgb.y)));
            image[ind+2] = min(255, int(round(rgb.z)));
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "You have to provide a scene file.\n";
        return 1;
    }

    Scene scene;
    scene.loadFromXml(argv[1]);
    for (const Mesh& mesh: scene.meshes) {
        for (const Face& face: mesh.faces) {
            scene.triangles.emplace_back(mesh.material_id, face);
        }
    }
    scene.meshes.clear();
    vector<Vec3f>& vertices = scene.vertex_data;
    for (Triangle& triangle: scene.triangles) {
        Vec3f& v0 = vertices[triangle.indices.v0_id - 1];
        Vec3f& v1 = vertices[triangle.indices.v1_id - 1];
        Vec3f& v2 = vertices[triangle.indices.v2_id - 1];
        Vec3f A = v1 - v0;
        Vec3f B = v2 - v0;
        triangle.normal = A.cross(B).normalize();
    }

    for (const Camera& camera: scene.cameras) {
        int width = camera.image_width;
        int height = camera.image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        Vec3f u = camera.up.cross(-camera.gaze).normalize();
        Vec3f middle = camera.position + (camera.gaze * camera.near_distance); // m
        Vec3f topLeft = middle + (u * camera.near_plane.x) + (camera.up * camera.near_plane.w); // q
        float pixelWidth = (camera.near_plane.y - camera.near_plane.x) / camera.image_width;
        float pixelHeight = (camera.near_plane.w - camera.near_plane.z) / camera.image_height;
        vector<thread> threads(4);
        int leap = height / 4;
        for (int i = 0; i < 4; ++i) {
            int l_height = i * leap;
            int r_height = (i == 3) ? height : (i + 1) * leap;
            threads[i] = thread(work, l_height, r_height, width, camera, topLeft, u, pixelWidth, pixelHeight, scene, image);
        }
        for (auto& t: threads) {
            t.join();
        }
        write_ppm(camera.image_name.c_str(), image, width, height);
        delete[] image;
    }
    return 0;
}
