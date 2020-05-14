//
// Created by Weugene on 2019-09-05.
//

#include <OpenGL/gl.h>

#include <OpenGl/glu.h>

#include <GLUT/glut.h>
//#include <gl/gl.h>
//#include <gl/glu.h>
void Render()
{
    //clear color and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();//load identity matrix

    glTranslatef(0.0f,0.0f,-4.0f);//move forward 4 units

    glColor3f(0.0f,0.0f,1.0f); //blue color

    glBegin(GL_TRIANGLES);//start drawing triangles
    glVertex3f(-1.0f,-0.25f,0.0f);//triangle one first vertex
    glVertex3f(-0.5f,-0.25f,0.0f);//triangle one second vertex
    glVertex3f(-0.75f,0.25f,0.0f);//triangle one third vertex
    //drawing a new triangle
    glVertex3f(0.5f,-0.25f,0.0f);//triangle two first vertex
    glVertex3f(1.0f,-0.25f,0.0f);//triangle two second vertex
    glVertex3f(0.75f,0.25f,0.0f);//triangle two third vertex
    glEnd();//end drawing of triangles
}

int main(){
    printf("Hello world!");
//    Render();
    return  0;
}