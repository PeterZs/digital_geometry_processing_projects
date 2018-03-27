#version 330

layout (location = 0) out vec4 FragColor;

uniform vec3 bot_left;
uniform vec3 top_right;

void main()
{
    float delta = 0.5;
    // selection is currently only 1px thick
    if (gl_FragCoord.x - bot_left.x     > delta &&
        top_right.x    - gl_FragCoord.x > delta &&
        gl_FragCoord.y - bot_left.y     > delta &&
        top_right.y    - gl_FragCoord.y > delta)
        FragColor = vec4(1.0, 0.0, 0.0, 0.1); // fill inner selection with transparent red
    else
        FragColor = vec4(1.0, 0.0, 0.0, 1.0); // draw selection in red
}