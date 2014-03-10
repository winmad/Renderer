uniform mat4 gl_ModelViewMatrix;
uniform mat4 gl_ProjectionMatrix;

void main()
{
    gl_Position =gl_ProjectionMatrix * gl_ModelViewMatrix * gl_Vertex;
}