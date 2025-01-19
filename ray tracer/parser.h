#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <cmath>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        float x, y, z;

        Vec3f(): x(0), y(0), z(0) {}

        Vec3f(float x, float y, float z) {
            this->x = x;
            this->y = y;
            this->z = z;
        }

        Vec3f operator*(float c) const {
            return Vec3f(c * x, c * y, c * z);
        }

        Vec3f operator+(Vec3f rhs) const {
            return Vec3f(x + rhs.x, y + rhs.y, z + rhs.z);
        }

        Vec3f operator+=(Vec3f rhs) {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;

            return *this;
        }

        Vec3f operator-() const {
            return Vec3f(-x, -y, -z);
        }

        Vec3f operator-(Vec3f rhs) const {
            return Vec3f(x - rhs.x, y - rhs.y, z - rhs.z);
        }

        float operator*(Vec3f rhs) const {
            return x * rhs.x + y * rhs.y + z * rhs.z;
        }

        Vec3f operator*=(float c) {
            x *= c;
            y *= c;
            z *= c;

            return *this;
        }

        Vec3f operator/(float size) const {
            if (size == 0) {
                return *this;
            }

            return Vec3f(x / size, y / size, z / size);
        }

        float square() const {
            return x * x + y * y + z * z;
        }

        Vec3f normalize() {
            float len = sqrt(square());
            if (len == 0) {
                return *this;
            }
            
            return *this / len;
        }

        float length() {
            return sqrt(square());
        }

        Vec3f cross(Vec3f rhs) const {
            return Vec3f(y * rhs.z - z * rhs.y, 
                        z * rhs.x - x * rhs.z, 
                        x * rhs.y - y * rhs.x);
        }
    };

    struct Vec3i
    {
        int x, y, z;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f operator-(Vec3f rhs) const {
            return Vec3f(position.x - rhs.x, position.y - rhs.y, position.z - rhs.z);
        } 

        Vec3f position;
        Vec3f intensity; // TODO
    };

    struct Material
    {
        bool is_mirror;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
    };

    struct Triangle
    {
        Triangle() {}
        Triangle(int material_id, Face indices) {
            this->material_id = material_id;
            this->indices = indices;
        }
        Vec3f normal;
        int material_id;
        Face indices;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
    };

    struct Scene
    {

        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        //Functions
        void loadFromXml(const std::string &filepath);
    };
}

#endif
