#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp> // GL Math library header
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp> 
#include <ft2build.h>
#include FT_FREETYPE_H

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;

GLuint gProgram[3];
int gWidth = 600, gHeight = 1000;
GLuint gVertexAttribBuffer, gTextVBO, gIndexBuffer;
GLuint gTex2D;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
int gTriangleIndexDataSizeInBytes, gLineIndexDataSizeInBytes;

GLint modelingMatrixLoc[2];
GLint viewingMatrixLoc[2];
GLint projectionMatrixLoc[2];
GLint eyePosLoc[2];
GLint lightPosLoc[2];
GLint kdLoc[2];

glm::mat4 projectionMatrix;
glm::mat4 viewingMatrix;
glm::mat4 modelingMatrix = glm::translate(glm::mat4(1.f), glm::vec3(-4.5f, -7.5f, -4.5f));
glm::vec3 eyePos = glm::vec3(0, 0, 26);
glm::vec3 lightPos = glm::vec3(0, 0, 7);

glm::vec3 kdGround(0.334, 0.288, 0.635); // this is the ground color in the demo
glm::vec3 kdCubes(0.86, 0.11, 0.31);

int activeProgramIndex = 0;

std::string facings[4] = { "Front", "Right", "Back", "Left" };
int facing_id = 0;
int points = 0;
std::string pressed_key = "";
double key_press_timestamp = -99;
int game_area[9][15][9];
double last_timestamp = glfwGetTime();
double time_between_fall = 1;
double current_angle = 0;
double destination_angle = 0;
bool game_over = false;

struct Point {
    int x, y, z;
    Point(int x, int y, int z): x(x), y(y), z(z) {}
    Point(): x(4), y(13), z(4) {}
};

Point active_cube_center;

void new_active_cube() {
    Point test_cube = Point();
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                int x = test_cube.x + i;
                int y = test_cube.y + j;
                int z = test_cube.z + k;
                if (game_area[x][y][z]) {
                    game_over = true;
                    return;
                }
            }
        }
    }
    active_cube_center = Point();

    // std::cout << active_cube_center.x << " " << active_cube_center.y << " " << active_cube_center.z << "\n";
}

bool check_collision_in_next_time_instant() {
    // hit ground
    if (active_cube_center.y - 1 == 0) {
        return true;
    }
    // hit another piece
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            int x_under = active_cube_center.x + i;
            int y_under = active_cube_center.y - 2;
            int z_under = active_cube_center.z + j;
            if (game_area[x_under][y_under][z_under]) {
                return true;
            }
        }
    }
    return false;
}

void shift_blocks_down(int cleared_level) {
    for (int j = cleared_level + 1; j < 15; ++j) {
        for (int i = 0; i < 9; ++i) {
            for (int k = 0; k < 9; ++k) {
                game_area[i][j - 1][k] = game_area[i][j][k];
                if (j == 14) {
                    game_area[i][j][k] = 0;
                }
            }
        }
    }
}

void clear_level() {
    for (int j = 0; j < 15; ++j) {
        bool all_ones = true;
        for (int i = 0; i < 9; ++i) {
            for (int k = 0; k < 9; ++k) {
                if (!game_area[i][j][k]) {
                    all_ones = false;
                    break;
                }
            }
            if (!all_ones) {
                break;
            }
        }
        if (all_ones) {
            points += 81;
            shift_blocks_down(j);
        }
    }
}

void mark_game_area() {
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                int x = active_cube_center.x + i;
                int y = active_cube_center.y + j;
                int z = active_cube_center.z + k;
                game_area[x][y][z] = 1;
            }
        }
    }
}

void go_down() {
    if (check_collision_in_next_time_instant()) {
        mark_game_area();
        new_active_cube();
        return;
    }
    active_cube_center.y -= 1;
}

void go_left() {
    if (facing_id == 0) {
        if (active_cube_center.x > 1) {
            int left_center = active_cube_center.x - 2;
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    if (game_area[left_center][active_cube_center.y + j][active_cube_center.z + k]) {
                        return;
                    }
                }
            }
            active_cube_center.x -= 1;
        }
    } else if (facing_id == 1) {
        if (active_cube_center.z < 7) {
            int left_center = active_cube_center.z + 2;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if (game_area[active_cube_center.x + i][active_cube_center.y + j][left_center]) {
                        return;
                    }
                }
            }
            active_cube_center.z += 1;
        }
    } else if (facing_id == 2) {
        if (active_cube_center.x < 7) {
            int left_center = active_cube_center.x + 2;
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    if (game_area[left_center][active_cube_center.y + j][active_cube_center.z + k]) {
                        return;
                    }
                }
            }
            active_cube_center.x += 1;
        }
    } else {
        if (active_cube_center.z > 1) {
            int left_center = active_cube_center.z - 2;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if (game_area[active_cube_center.x + i][active_cube_center.y + j][left_center]) {
                        return;
                    }
                }
            }
            active_cube_center.z -= 1;
        }
    }
}

void go_right() {
    if (facing_id == 0) {
        if (active_cube_center.x < 7) {
            int right_center = active_cube_center.x + 2;
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    if (game_area[right_center][active_cube_center.y + j][active_cube_center.z + k]) {
                        return;
                    }
                }
            }
            active_cube_center.x += 1;
        }
    } else if (facing_id == 1) {
        if (active_cube_center.z > 1) {
            int right_center = active_cube_center.z - 2;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if (game_area[active_cube_center.x + i][active_cube_center.y + j][right_center]) {
                        return;
                    }
                }
            }
            active_cube_center.z -= 1;
        }
    } else if (facing_id == 2) {
        if (active_cube_center.x > 1) {
            int right_center = active_cube_center.x - 2;
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    if (game_area[right_center][active_cube_center.y + j][active_cube_center.z + k]) {
                        return;
                    }
                }
            }
            active_cube_center.x -= 1;
        }
    } else {
        if (active_cube_center.z < 7) {
            int right_center = active_cube_center.z + 2;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if (game_area[active_cube_center.x + i][active_cube_center.y + j][right_center]) {
                        return;
                    }
                }
            }
            active_cube_center.z += 1;
        }
    }
}


// Holds all state information relevant to a character as loaded using FreeType
struct Character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

// For reading GLSL files
bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

GLuint createVS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

	return vs;
}

GLuint createFS(const char* shaderName)
{
    string shaderSource;

    string filename(shaderName);
    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

	return fs;
}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[2]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", 0, &face))
    //if (FT_New_Face(ft, "/usr/share/fonts/truetype/gentium-basic/GenBkBasR.ttf", 0, &face)) // you can use different fonts
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
                );
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            (GLuint) face->glyph->advance.x
        };
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void initShaders()
{
	// Create the programs

    gProgram[0] = glCreateProgram();
	gProgram[1] = glCreateProgram();
	gProgram[2] = glCreateProgram();

	// Create the shaders for both programs

    GLuint vs1 = createVS("vert.glsl"); // for cube shading
    GLuint fs1 = createFS("frag.glsl");

	GLuint vs2 = createVS("vert2.glsl"); // for border shading
	GLuint fs2 = createFS("frag2.glsl");

	GLuint vs3 = createVS("vert_text.glsl");  // for text shading
	GLuint fs3 = createFS("frag_text.glsl");

	// Attach the shaders to the programs

	glAttachShader(gProgram[0], vs1);
	glAttachShader(gProgram[0], fs1);

	glAttachShader(gProgram[1], vs2);
	glAttachShader(gProgram[1], fs2);

	glAttachShader(gProgram[2], vs3);
	glAttachShader(gProgram[2], fs3);

	// Link the programs

    for (int i = 0; i < 3; ++i)
    {
        glLinkProgram(gProgram[i]);
        GLint status;
        glGetProgramiv(gProgram[i], GL_LINK_STATUS, &status);

        if (status != GL_TRUE)
        {
            cout << "Program link failed: " << i << endl;
            exit(-1);
        }
    }


	// Get the locations of the uniform variables from both programs

	for (int i = 0; i < 2; ++i)
	{
		modelingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "modelingMatrix");
		viewingMatrixLoc[i] = glGetUniformLocation(gProgram[i], "viewingMatrix");
		projectionMatrixLoc[i] = glGetUniformLocation(gProgram[i], "projectionMatrix");
		eyePosLoc[i] = glGetUniformLocation(gProgram[i], "eyePos");
		lightPosLoc[i] = glGetUniformLocation(gProgram[i], "lightPos");
		kdLoc[i] = glGetUniformLocation(gProgram[i], "kd");

        glUseProgram(gProgram[i]);
        glUniformMatrix4fv(modelingMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(modelingMatrix));
        glUniform3fv(eyePosLoc[i], 1, glm::value_ptr(eyePos));
        glUniform3fv(lightPosLoc[i], 1, glm::value_ptr(lightPos));
        glUniform3fv(kdLoc[i], 1, glm::value_ptr(kdCubes));
	}
}

// VBO setup for drawing a cube and its borders
void initVBO()
{
    GLuint vao;
    glGenVertexArrays(1, &vao);
    assert(vao > 0);
    glBindVertexArray(vao);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	assert(glGetError() == GL_NONE);

	glGenBuffers(1, &gVertexAttribBuffer);
	glGenBuffers(1, &gIndexBuffer);

	assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    GLuint indices[] = {
        0, 1, 2, // front
        3, 0, 2, // front
        4, 7, 6, // back
        5, 4, 6, // back
        0, 3, 4, // left
        3, 7, 4, // left
        2, 1, 5, // right
        6, 2, 5, // right
        3, 2, 7, // top
        2, 6, 7, // top
        0, 4, 1, // bottom
        4, 5, 1  // bottom
    };

    GLuint indicesLines[] = {
        7, 3, 2, 6, // top
        4, 5, 1, 0, // bottom
        2, 1, 5, 6, // right
        5, 4, 7, 6, // back
        0, 1, 2, 3, // front
        0, 3, 7, 4, // left
    };

    GLfloat vertexPos[] = {
        0, 0, 1, // 0: bottom-left-front
        1, 0, 1, // 1: bottom-right-front
        1, 1, 1, // 2: top-right-front
        0, 1, 1, // 3: top-left-front
        0, 0, 0, // 0: bottom-left-back
        1, 0, 0, // 1: bottom-right-back
        1, 1, 0, // 2: top-right-back
        0, 1, 0, // 3: top-left-back
    };

    GLfloat vertexNor[] = {
         1.0,  1.0,  1.0, // 0: unused
         0.0, -1.0,  0.0, // 1: bottom
         0.0,  0.0,  1.0, // 2: front
         1.0,  1.0,  1.0, // 3: unused
        -1.0,  0.0,  0.0, // 4: left
         1.0,  0.0,  0.0, // 5: right
         0.0,  0.0, -1.0, // 6: back 
         0.0,  1.0,  0.0, // 7: top
    };

	gVertexDataSizeInBytes = sizeof(vertexPos);
	gNormalDataSizeInBytes = sizeof(vertexNor);
    gTriangleIndexDataSizeInBytes = sizeof(indices);
    gLineIndexDataSizeInBytes = sizeof(indicesLines);
    int allIndexSize = gTriangleIndexDataSizeInBytes + gLineIndexDataSizeInBytes;

	glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexPos);
	glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, vertexNor);

	glBufferData(GL_ELEMENT_ARRAY_BUFFER, allIndexSize, 0, GL_STATIC_DRAW);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, gTriangleIndexDataSizeInBytes, indices);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, gTriangleIndexDataSizeInBytes, gLineIndexDataSizeInBytes, indicesLines);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));
}

void init() 
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    // polygon offset is used to prevent z-fighting between the cube and its borders
    glPolygonOffset(0.5, 0.5);
    glEnable(GL_POLYGON_OFFSET_FILL);

    initShaders();
    initVBO();
    initFonts(gWidth, gHeight);
}

void drawCube()
{
	glUseProgram(gProgram[0]);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

void drawCubeEdges()
{
    glLineWidth(3);

	glUseProgram(gProgram[1]);

    for (int i = 0; i < 6; ++i)
    {
	    glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_INT, BUFFER_OFFSET(gTriangleIndexDataSizeInBytes + i * 4 * sizeof(GLuint)));
    }
}

void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    // Activate corresponding render state	
    glUseProgram(gProgram[2]);
    glUniform3f(glGetUniformLocation(gProgram[2], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },            
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }           
        };

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}

void display_active_cube() {
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                int x = active_cube_center.x + i;
                int y = active_cube_center.y + j;
                int z = active_cube_center.z + k;
                glm::mat4 cur_modeling_matrix = modelingMatrix * glm::translate(glm::mat4(1.f), glm::vec3(x, y, z));
                for (int l = 0; l < 2; ++l) {
                    glUseProgram(gProgram[l]);
                    glUniformMatrix4fv(modelingMatrixLoc[l], 1, GL_FALSE, glm::value_ptr(cur_modeling_matrix));
                    glUniform3fv(kdLoc[l], 1, glm::value_ptr(kdCubes));
                }
                drawCube();
                drawCubeEdges();
            }
        }
    }
}

void display_game_area() {
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 15; ++j) {
            for (int k = 0; k < 9; ++k) {
                if (game_area[i][j][k]) {
                    glm::mat4 cur_modeling_matrix = modelingMatrix * glm::translate(glm::mat4(1.f), glm::vec3(i, j, k));
                    for (int l = 0; l < 2; ++l) {
                        glUseProgram(gProgram[l]);
                        glUniformMatrix4fv(modelingMatrixLoc[l], 1, GL_FALSE, glm::value_ptr(cur_modeling_matrix));
                        glUniform3fv(kdLoc[l], 1, glm::value_ptr(kdCubes));
                    }
                    drawCube();
                    drawCubeEdges();
                }
            }
        }
    }
}

void display_ground() {
    for (int i = 0; i < 9; ++i) {
        for (int k = 0; k < 9; ++k) {
            glm::mat4 cur_modeling_matrix = modelingMatrix * glm::translate(glm::mat4(1.f), glm::vec3(i, -0.5f, k)) * glm::scale(glm::mat4(1.f), glm::vec3(1.f, 0.5f, 1.f));
            for (int l = 0; l < 2; ++l) {
                glUseProgram(gProgram[l]);
                glUniformMatrix4fv(modelingMatrixLoc[l], 1, GL_FALSE, glm::value_ptr(cur_modeling_matrix));
                glUniform3fv(kdLoc[l], 1, glm::value_ptr(kdGround));
            }
            drawCube();
            drawCubeEdges();
        }
    }
}

void rotate_camera_and_light() {
    // always look toward (0, 0, 0)
    float cur_angle_rad = (current_angle / 180.) * M_PI;
    glm::vec3 cur_eye_pos = glm::rotate<float>(glm::mat4(1.f), cur_angle_rad, glm::vec3(0.f, 1.f, 0.f)) * glm::vec4(eyePos, 1.f);
	viewingMatrix = glm::lookAt(cur_eye_pos, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));
    glm::vec3 cur_light_pos = glm::rotate<float>(glm::mat4(1.f), cur_angle_rad, glm::vec3(0.f, 1.f, 0.f)) * glm::vec4(lightPos, 1.f);

    for (int i = 0; i < 2; ++i)
    {
        glUseProgram(gProgram[i]);
        glUniformMatrix4fv(viewingMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
        glUniform3fv(lightPosLoc[i], 1, glm::value_ptr(cur_light_pos));
    }
}

void display()
{
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    // display all things
    display_game_area();
    display_active_cube();
    display_ground();

    // rotate camera and light
    if (current_angle < destination_angle) {
        current_angle += 2;
    }
    if (current_angle > destination_angle) {
        current_angle -= 2;
    }
    rotate_camera_and_light();

    // go down by one square and check if level is cleared
    if (!game_over) {
        if (time_between_fall != -1 && glfwGetTime() - last_timestamp >= time_between_fall) {
            last_timestamp = glfwGetTime();
            go_down();
        }
        clear_level();
    }

    // display text
    renderText(facings[facing_id], 25, gHeight - 40, 0.75, glm::vec3(1, 1, 0));
    if (glfwGetTime() - key_press_timestamp <= 0.6) {
        renderText(pressed_key, 25, gHeight - 90, 0.75, glm::vec3(1, 0, 0));
    }
    if (game_over) {
        renderText("Game Over!", gWidth / 2 - 220, gHeight / 2, 1.75, glm::vec3(1, 1, 0));
    }
    renderText("Points: " + std::to_string(points), gWidth - 200, gHeight - 40, 0.75, glm::vec3(1, 1, 0));

    assert(glGetError() == GL_NO_ERROR);
}

void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);

	// Use perspective projection

	float fovyRad = (float) (45.0 / 180.0) * M_PI;
	projectionMatrix = glm::perspective(fovyRad, gWidth / (float) gHeight, 1.0f, 100.0f);

    // always look toward (0, 0, 0)
	viewingMatrix = glm::lookAt(eyePos, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));

    for (int i = 0; i < 2; ++i)
    {
        glUseProgram(gProgram[i]);
        glUniformMatrix4fv(projectionMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(projectionMatrix));
        glUniformMatrix4fv(viewingMatrixLoc[i], 1, GL_FALSE, glm::value_ptr(viewingMatrix));
    }
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action != GLFW_PRESS) {
        return;
    }
    if (key == GLFW_KEY_Q || key == GLFW_KEY_ESCAPE) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
        return;
    }
    if (game_over) {
        return;
    }
    switch (key) {
        case GLFW_KEY_A:
            go_left();
            pressed_key = "A";
            key_press_timestamp = glfwGetTime();
            break;
        case GLFW_KEY_D:
            go_right();
            pressed_key = "D";
            key_press_timestamp = glfwGetTime();
            break;
        case GLFW_KEY_S:
            if (time_between_fall == -1) {
                time_between_fall = 1;
            }
            if (time_between_fall > 0.4) {
                time_between_fall -= 0.2;
            }
            pressed_key = "S";
            key_press_timestamp = glfwGetTime();
            break;
        case GLFW_KEY_W:
            if (time_between_fall == 1) {
                time_between_fall = -1; // stop 
            } else if (time_between_fall != -1) {
                time_between_fall += 0.2;
            }
            pressed_key = "W";
            key_press_timestamp = glfwGetTime();
            break;
        case GLFW_KEY_H:
            facing_id = (facing_id + 3) % 4;
            destination_angle -= 90;
            pressed_key = "H";
            key_press_timestamp = glfwGetTime();
            break;
        case GLFW_KEY_K:
            facing_id = (facing_id + 1) % 4;
            destination_angle += 90;
            pressed_key = "K";
            key_press_timestamp = glfwGetTime();
            break;
    }
}

void mainLoop(GLFWwindow* window)
{
    new_active_cube();
    while (!glfwWindowShouldClose(window))
    {
        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char** argv) {
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(gWidth, gHeight, "tetrisGL", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, gWidth, gHeight); // need to call this once ourselves
    mainLoop(window); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}
