#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
    this->depth_buffer = std::vector<std::vector<double>>(camera->horRes, std::vector<double>(camera->verRes, 31.0));

	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}

Vec4 Scene::scale(Scaling& scaling, Vec4& v) {
    double sx = scaling.sx, sy = scaling.sy, sz = scaling.sz;
    double s[4][4] = {
            {sx,0,0,0},
            {0,sy,0,0},
            {0,0,sz,0},
            {0,0,0,1}
    };

    Matrix4 S = Matrix4(s);

    return multiplyMatrixWithVec4(S, v);
}

Vec4 Scene::rotate(Rotation& rotation, Vec4& v) {
    double theta = rotation.angle, ux = rotation.ux, uy = rotation.uy, uz = rotation.uz;
    theta *= M_PI / 180;

    Vec3 U(ux, uy, uz);
    Vec3 V;
    if (abs(ux) < abs(uy) && abs(ux) < abs(uy)) {
        V.x = 0;
        V.y = -uz;
        V.z = uy;
    } else if (abs(uy) < abs(uz)) {
        V.x = -uz;
        V.y = 0;
        V.z = ux;
    } else {
        V.x = -uy;
        V.y = ux;
        V.z = 0;
    }
    Vec3 W = crossProductVec3(U, V);
    U = normalizeVec3(U);
    V = normalizeVec3(V);
    W = normalizeVec3(W);

    double cameraTransformation[4][4] = {
            {U.x, U.y, U.z, 0},
            {V.x, V.y, V.z, 0},
            {W.x, W.y, W.z, 0},
            {0, 0, 0, 1},
    };

    double rotationTransformation[4][4] = {
            {1, 0, 0, 0},
            {0, cos(theta), -sin(theta), 0},
            {0, sin(theta), cos(theta)},
            {0, 0, 0, 1}
    };

    double cameraTransformationInverse[4][4] = {
            {U.x, V.x, W.x, 0},
            {U.y, V.y, W.y, 0},
            {U.z, V.z, W.z, 0},
            {0, 0, 0, 1},
    };

    Matrix4 CameraTransformation = Matrix4(cameraTransformation);
    Matrix4 RotationTransformation = Matrix4(rotationTransformation);
    Matrix4 CameraTransformationInverse = Matrix4(cameraTransformationInverse);

    Matrix4 CompositeTransformation = multiplyMatrixWithMatrix(
            CameraTransformationInverse,
            multiplyMatrixWithMatrix(RotationTransformation, CameraTransformation)
    );

    return multiplyMatrixWithVec4(CompositeTransformation, v);
}

Vec4 Scene::translate(Translation& translation, Vec4& v) {
    double tx = translation.tx, ty = translation.ty, tz = translation.tz;
    double t[4][4] = {
            {1,0,0,tx},
            {0,1,0,ty},
            {0,0,1,tz},
            {0,0,0,1}
    };

    Matrix4 T = Matrix4(t);
    return multiplyMatrixWithVec4(T, v);
}

Vec4 Scene::orthographicProjection(Camera* camera, Vec4 &vec) {
    double r = camera->right, l = camera->left, t = camera->top, b = camera->bottom, f = camera->far, n = camera->near;
    Vec3 u = camera->u, v = camera->v, w = camera->w, e = camera->position;

    double mcam[4][4] = {
            { u.x, u.y, u.z, -(u.x * e.x + u.y * e.y + u.z * e.z) },
            { v.x, v.y, v.z, -(v.x * e.x + v.y * e.y + v.z * e.z) },
            { w.x, w.y, w.z, -(w.x * e.x + w.y * e.y + w.z * e.z) },
            { 0, 0, 0, 1 }
    };
    Matrix4 Mcam(mcam);

    double morth[4][4] = {
        { 2.0 / (r - l), 0, 0, -(r + l) / (r - l) },
        { 0, 2.0 / (t - b), 0, -(t + b) / (t - b) },
        {0, 0, -2.0 / (f - n), -(f + n) / (f - n) },
        { 0, 0, 0, 1 }
    };
    Matrix4 Morth(morth);
    return multiplyMatrixWithVec4(Morth, multiplyMatrixWithVec4(Mcam, vec));
}

Vec4 Scene::perspectiveProjection(Camera* camera, Vec4 &vec) {
    double r = camera->right, l = camera->left, t = camera->top, b = camera->bottom, f = camera->far, n = camera->near;
    Vec3 u = camera->u, v = camera->v, w = camera->w, e = camera->position;

    double mcam[4][4] = {
            { u.x, u.y, u.z, -(u.x * e.x + u.y * e.y + u.z * e.z) },
            { v.x, v.y, v.z, -(v.x * e.x + v.y * e.y + v.z * e.z) },
            { w.x, w.y, w.z, -(w.x * e.x + w.y * e.y + w.z * e.z) },
            { 0, 0, 0, 1 }
    };
    Matrix4 Mcam(mcam);

    double mper[4][4] = {
            { 2 * n / (r - l), 0, (r + l) / (r - l), 0 },
            { 0, 2 * n / (t - b), (t + b) / (t - b), 0 },
            { 0, 0, -(f + n) / (f - n), -2 * f * n / (f - n) },
            { 0, 0, -1, 0 }
    };
    Matrix4 Mper(mper);
    Vec4 projected = multiplyMatrixWithVec4(Mper, multiplyMatrixWithVec4(Mcam, vec));
    // perspective divide
    projected.x /= projected.t;
    projected.y /= projected.t;
    projected.z /= projected.t;
    projected.t /= projected.t;
    return projected;
}

Vec4 Scene::viewportTransformation(Camera *camera, Vec4 &v) {
    double nx = camera->horRes, ny = camera->verRes;
    double mvp[4][4] = {
            { nx / 2, 0, 0, (nx - 1) / 2 },
            { 0, ny / 2, 0, (ny - 1) / 2 },
            { 0, 0, 0.5, 0.5 },
            { 0, 0, 0, 0 },
    };
    Matrix4 Mvp(mvp);
    return multiplyMatrixWithVec4(Mvp, v);
}

bool Scene::visible(double denom, double num, double& t_E, double& t_L) {
    if (denom > 0) {
        double t = num / denom;
        if (t > t_L) {
            return false;
        }

        if (t > t_E) {
            t_E = t;
        }
    }

    else if (denom < 0) {
        double t = num / denom;
        if (t < t_E) {
            return false;
        }

        if (t < t_L) {
            t_L = t;
        }
    }

    else if (num > 0) { // denom == 0
        return false;
    }

    return true;
}

bool Scene::clipLine(Vec4& v0, Vec4& v1, Color& c0, Color& c1) {
    double dx = v1.x - v0.x;
    double dy = v1.y - v0.y;
    double dz = v1.z - v0.z;
    double dcr = c1.r - c0.r;
    double dcg = c1.g - c0.g;
    double dcb = c1.b - c0.b;

    double dd = sqrt(pow(v1.x - v0.x, 2) + pow(v1.y - v0.y, 2) + pow(v1.z - v0.z, 2));

    double x_min = -1, y_min = -1, z_min = -1;
    double x_max = 1, y_max = 1, z_max = 1;

    double t_E = 0, t_L = 1;
    bool isVisible = false;

    double dist =  sqrt(dx * dx + dy * dy + dz * dz);

    // printf("%f %f %f\n", v0.x, v0.y, v0.z);
    // printf("%f %f %f\n", v1.x, v1.y, v1.z);

    if (visible(dx, x_min - v0.x, t_E, t_L)) {
        if (visible(-dx, v0.x - x_max, t_E, t_L)) {
            if (visible(dy, y_min - v0.y, t_E, t_L)) {
                if (visible(-dy, v0.y - y_max, t_E, t_L)) {
                    if (visible(dz, z_min - v0.z, t_E, t_L)) {
                        if (visible(-dz, v0.z - z_max, t_E, t_L)) {
                            isVisible = true;

                            if (t_L < 1) {
                                v1.x = v0.x + dx * t_L;
                                v1.y = v0.y + dy * t_L;
                                v1.z = v0.z + dz * t_L;

                                c1.r = c0.r + (dcr / dd) * dist * t_L;
                                c1.g = c0.g + (dcg / dd) * dist * t_L;
                                c1.b = c0.b + (dcb / dd) * dist * t_L;
                                // printf("1 %f %f %f %f\n", dx, dy, dz, t_L);
                                // printf("%f %f %f\n", v0.x, v0.y, v0.z);
                                // printf("%f %f %f\n", v1.x, v1.y, v1.z);
                            }

                            if (t_E > 0) {
                                v0.x += dx * t_E;
                                v0.y += dy * t_E;
                                v0.z += dz * t_E;

                                c0.r += (dcr / dd) * dist * t_E;
                                c1.g += (dcg / dd) * dist * t_E;
                                c1.b += (dcb / dd) * dist * t_E;

                                // printf("2 %f %f %f %f\n", dx, dy, dz, t_E);
                                // printf("%f %f %f\n", v0.x, v0.y, v0.z);
                                // printf("%f %f %f\n", v1.x, v1.y, v1.z);
                            }
                        }
                    }
                }
            }
        }
    }

    return isVisible;
}

void Scene::drawLine(Vec4& v0, Vec4& v1, Color& c0, Color& c1) {
    double m;
    if (v0.x == v1.x) {
        if (v0.y < v1.y) {
            m = FLT_MAX;
        }
        else {
            m = FLT_MIN;
        }
    }
    else {
        m = (v1.y - v0.y) / (v1.x - v0.x);
    }
    // printf("%f\n", m);
    if (v0.x > v1.x) { // v0'ı sola at
        swap(v0, v1);
        swap(c0, c1);
    }
    if (0.0 < m && m <= 1.0) {
        double y = v0.y;
        double d = (v0.y - v1.y) + 0.5 * (v1.x - v0.x);

        double increment = (v0.y - v1.y) + (v1.x - v0.x);
        double decrement = (v0.y - v1.y);

        for (int x = v0.x; x <= v1.x; x++) {
            Color c = interpolateColor(x, v0, v1, c0, c1, 'x');
            double depth = interpolateDepth(x, v0, v1, 'x');
            draw(x, y, depth, c);

            if (d < 0) { // NE
                y++;
                d += increment;
            }
            else { // E
                d += decrement;
            }
        }
    }

    else if (m > 1.0) {
        double x = v0.x;
        double d = (v0.x - v1.x) + 0.5 * (v1.y - v0.y);

        double increment = (v0.x - v1.x) + (v1.y - v0.y);
        double decrement = (v0.x - v1.x);

        for (int y = v0.y; y <= v1.y; y++) {
            Color c = interpolateColor(y, v0, v1, c0, c1, 'y');
            double depth = interpolateDepth(y, v0, v1, 'y');
            draw(x, y, depth, c);

            if (d < 0) {
                x++;
                d += increment;
            }
            else {
                d += decrement;
            }
        }
    }

    else if (-1.0 <= m && m <= 0) {
        double y = v1.y;
        double d = -(v0.y - v1.y) + 0.5 * (v1.x - v0.x);

        double increment = -(v0.y - v1.y) + (v1.x - v0.x);
        double decrement = -(v0.y - v1.y);

        for (int x = v1.x; x >= v0.x; x--) {
            Color c = interpolateColor(x, v0, v1, c0, c1, 'x');
            double depth = interpolateDepth(x, v0, v1, 'x');
            draw(x, y, depth, c);

            if (d < 0) { // NW
                y++;
                d += increment;
            }
            else { // W
                d += decrement;
            }
        }
    }

    else if (m < -1.0) {
        double x = v0.x;
        double d = (v0.x - v1.x) - 0.5 * (v1.y - v0.y);

        double increment = (v0.x - v1.x) - (v1.y - v0.y);
        double decrement = (v0.x - v1.x);

        for (int y = v0.y; y >= v1.y; y--) {
            Color c = interpolateColor(y, v0, v1, c0, c1, 'y');
            double depth = interpolateDepth(y, v0, v1, 'y');
            draw(x, y, depth, c);

            if (d < 0) {
                x++;
                d += increment;
            }
            else {
                d += decrement;
            }
        }
    }
}

double Scene::f(double x, double y, double x0, double y0, double x1, double y1) {
    return x * (y0 - y1) + y * (x1 - x0) + x0 * y1 - y0 * x1;
}

void Scene::triangleRasterization(std::vector<Vec4>& triangle, Camera* camera) {
    Vec4 v0 = triangle[0], v1 = triangle[1], v2 = triangle[2];
    int x_min = min(min(v0.x, v1.x), v2.x);
    x_min = max(0, x_min);
    int x_max = max(max(v0.x, v1.x), v2.x);
    x_max = min(camera->horRes - 1, x_max);
    int y_min = min(min(v0.y, v1.y), v2.y);
    y_min = max(0, y_min);
    int y_max = max(max(v0.y, v1.y), v2.y);
    y_max = min(camera->verRes - 1, y_max);
    for (int y = y_min; y <= y_max; ++y) {
        for (int x = x_min; x <= x_max; ++x) {
            double alpha = f(x, y, v1.x, v1.y, v2.x, v2.y) / f(v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
            double beta = f(x, y, v2.x, v2.y, v0.x, v0.y) / f(v1.x, v1.y, v2.x, v2.y, v0.x, v0.y);
            double gamma = f(x, y, v0.x, v0.y, v1.x, v1.y) / f(v2.x, v2.y, v0.x, v0.y, v1.x, v1.y);
            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                double z_val = alpha * v0.z + beta * v1.z + gamma * v2.z;
                if (0 <= z_val && z_val <= 1 && z_val < this->depth_buffer[x][y]) {
                    this->depth_buffer[x][y] = z_val;
                    Color A = this->multiplyColor(*this->colorsOfVertices[v0.colorId - 1], alpha);
                    Color B = this->multiplyColor(*this->colorsOfVertices[v1.colorId - 1], beta);
                    Color C = this->multiplyColor(*this->colorsOfVertices[v2.colorId - 1], gamma);
                    this->image[x][y].r = A.r + B.r + C.r;
                    this->image[x][y].g = A.g + B.g + C.g;
                    this->image[x][y].b = A.b + B.b + C.b;
                }
            }
        }
    }
}

double Scene::interpolateDepth(int p, Vec4& v0, Vec4& v1, char type) {
    double alpha;
    if (type == 'x') {
        alpha = (p - v0.x) / (v1.x - v0.x);
    }
    else if (type == 'y') {
        alpha = (p - v0.y) / (v1.y - v0.y);
    }

    double term1 = v0.z * (1 - alpha);
    double term2 = v1.z * alpha;

    return term1 + term2;
}

Color Scene::interpolateColor(int p, Vec4& v0, Vec4& v1, Color& c0, Color& c1, char type) {
    double alpha;
    if (type == 'x') {
        alpha = (p - v0.x) / (v1.x - v0.x);
    } else if (type == 'y') {
        alpha = (p - v0.y) / (v1.y - v0.y);
    }

    Color term1 = this->multiplyColor(c0, 1 - alpha);
    Color term2 = this->multiplyColor(c1, alpha);

    return Color(term1.r + term2.r, term1.g + term2.g, term1.b + term2.b);
}

Color Scene::multiplyColor(Color color, double val) {
    color.r *= val;
    color.g *= val;
    color.b *= val;
    return color;
}

void Scene::draw(int x, int y, double depth, Color& c) {
    if (x < 0 || y < 0 || x >= this->image.size() || y >= this->image[0].size()) {
        return;
    }

    if (0 <= depth && depth <= 1 && depth < this->depth_buffer[x][y]) {
        this->depth_buffer[x][y] = depth;
        this->image[x][y].r = c.r;
        this->image[x][y].g = c.g;
        this->image[x][y].b = c.b;
    }
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera) {
    for (const auto& mesh: meshes) {
        int numberOfTransformations = mesh->numberOfTransformations;

        // modeling transformations
        for (const auto& triangle: mesh->triangles) {
            std::vector<Vec4> transformed_triangle;
            for (const auto& vertexId: triangle.vertexIds) {
                Vec3* vertex = this->vertices[vertexId - 1];
                Vec4 homogeneous_vertex(vertex->x, vertex->y, vertex->z, 1, vertex->colorId);

                for (int j = 0; j < numberOfTransformations; ++j) {
                    char transformationType = mesh->transformationTypes[j];
                    int transformationId = mesh->transformationIds[j] - 1;
                    if (transformationType == 'r') {
                        homogeneous_vertex = this->rotate(*this->rotations[transformationId], homogeneous_vertex);
                    } else if (transformationType == 't') {
                        homogeneous_vertex = this->translate(*this->translations[transformationId], homogeneous_vertex);
                    } else if (transformationType == 's') {
                        homogeneous_vertex = this->scale(*this->scalings[transformationId], homogeneous_vertex);
                    } else {
                        printf("Invalid transformation type.\n");
                        exit(1);
                    }
                }
                transformed_triangle.push_back(homogeneous_vertex);
            }

            // viewing transformations
            for (auto& vertex: transformed_triangle) {
                if (camera->projectionType == ORTOGRAPHIC_PROJECTION) {
                    vertex = this->orthographicProjection(camera, vertex);
                } else if (camera->projectionType == PERSPECTIVE_PROJECTION) {
                    vertex = this->perspectiveProjection(camera, vertex);
                } else {
                    printf("Invalid projection type.\n");
                    exit(1);
                }
            }

            if (this->cullingEnabled) {
                Vec4 &v0 = transformed_triangle[0], &v1 = transformed_triangle[1], &v2 = transformed_triangle[2];
                Vec3 A = Vec3(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
                Vec3 B = Vec3(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
                Vec3 N = crossProductVec3(A, B);

                Vec3 G = Vec3((v0.x + v1.x + v2.x) / 3, (v0.y + v1.y + v2.y) / 3, (v0.z + v1.z + v2.z) / 3);

                // printf("%f %s\n", dotProductVec3(N, G), camera->outputFilename.c_str());
                // printf("%d %d %d\n", triangle.vertexIds[0], triangle.vertexIds[1], triangle.vertexIds[2]);
                if (dotProductVec3(N, G) < 0) {
                    continue;
                }
            }

            // clipping - viewport transformation - rasterization
            if (mesh->type == WIREFRAME_MESH) {
                Vec4 &v0 = transformed_triangle[0], &v1 = transformed_triangle[1], &v2 = transformed_triangle[2];
                Vec4 new_v0, new_v1, new_v2;
                Color c0, c1, c2;

                c0 = *this->colorsOfVertices[v0.colorId - 1];
                c1 = *this->colorsOfVertices[v1.colorId - 1];
                new_v0 = v0;
                new_v1 = v1;
                if (clipLine(new_v0, new_v1, c0, c1)) {
                    new_v0 = viewportTransformation(camera, new_v0);
                    new_v1 = viewportTransformation(camera, new_v1);
                    // printf("%f %f %f\n", new_v0.x, new_v0.y, new_v0.z);
                    // printf("%f %f %f\n", new_v1.x, new_v1.y, new_v1.z);
                    drawLine(new_v0, new_v1, c0, c1);
                }


                c1 = *this->colorsOfVertices[v1.colorId - 1];
                c2 = *this->colorsOfVertices[v2.colorId - 1];
                new_v1 = v1;
                new_v2 = v2;
                if (clipLine(new_v1, new_v2, c1, c2)) {
                    new_v1 = viewportTransformation(camera, new_v1);
                    new_v2 = viewportTransformation(camera, new_v2);
                    drawLine(new_v1, new_v2, c1, c2);
                }


                c2 = *this->colorsOfVertices[v2.colorId - 1];
                c0 = *this->colorsOfVertices[v0.colorId - 1];
                new_v2 = v2;
                new_v0 = v0;
                if (clipLine(new_v2, new_v0, c2, c0)) {
                    new_v2 = viewportTransformation(camera, new_v2);
                    new_v0 = viewportTransformation(camera, new_v0);
                    drawLine(new_v2, new_v0, c2, c0);
                }

            } else if (mesh->type == SOLID_MESH) {
                for (auto& vertex: transformed_triangle) {
                    vertex = viewportTransformation(camera, vertex);
                }
                triangleRasterization(transformed_triangle, camera);
            } else {
                printf("Invalid mesh type.\n");
                exit(1);
            }
        }
    }
}
