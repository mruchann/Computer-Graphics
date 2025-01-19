#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

    std::vector<std::vector<double>> depth_buffer;
	std::vector<std::vector<Color>> image;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);
	void forwardRenderingPipeline(Camera *camera);

    Vec4 scale(Scaling& scaling, Vec4& v);
    Vec4 rotate(Rotation& rotation, Vec4& v);
    Vec4 translate(Translation& translation, Vec4& v);
    Vec4 perspectiveProjection(Camera* camera, Vec4& v);
    Vec4 orthographicProjection(Camera* camera, Vec4& v);
    Vec4 viewportTransformation(Camera* camera, Vec4& v);

    bool visible(double denom, double num, double& t_E, double& t_L);
    bool clipLine(Vec4& v0, Vec4& v1, Color& c0, Color& c1);
    void drawLine(Vec4& v0, Vec4& v1, Color& c0, Color& c1);

    double f(double x, double y, double x0, double y0, double x1, double y1);
    void triangleRasterization(std::vector<Vec4>& triangle, Camera* camera);

    double interpolateDepth(int p, Vec4& v0, Vec4& v1, char type);
    Color interpolateColor(int p, Vec4& v0, Vec4& v1, Color& c0, Color& c1, char type);
    Color multiplyColor(Color color, double val);

    void draw(int x, int y, double depth, Color& c);
};

#endif
