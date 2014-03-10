#extension GL_EXT_gpu_shader4 : enable

uniform int nTriangles;
uniform int nQueries;

uniform sampler2D vertexPositionTex;
uniform sampler2D triangleVertexIndicesTex;

uniform sampler2D queryRayOriginTex;
uniform sampler2D queryRayDirectionTex;

float intersectTriangle(vec3 v1, vec3 v2, vec3 v3, vec3 origin, vec3 direction);

vec4 fetch(sampler2D tex, int index);

void main()
{
	int queryIndex = int(gl_FragCoord.y) * textureSize(queryRayOriginTex, 0).x + int(gl_FragCoord.x);
	if(queryIndex >= nQueries)
	{
		gl_FragColor = vec4(nQueries);
		return;
	}
	vec3 origin = fetch(queryRayOriginTex, queryIndex).xyz;
	vec3 direction = fetch(queryRayDirectionTex, queryIndex).xyz;
	vec4 dist_tid = vec4(-1, -1, 0, 0);
	
	for(int ti=0; ti<nTriangles; ti++)
	{
		vec3 indices = fetch(triangleVertexIndicesTex, ti).xyz;
		vec3 vertices[3], lines[3], r2v, n;
		for(int vi=0; vi<3; vi++)
			vertices[vi] = fetch(vertexPositionTex, int(indices[vi])).xyz;
		float dist = intersectTriangle(vertices[0], vertices[1], vertices[2], origin, direction);
		if(dist >= 0 && (dist_tid.x < 0 || dist_tid.x > dist))
		{
			dist_tid = vec4(dist, ti, 0, 0);
		}
	}
	gl_FragColor = vec4(dist_tid.x-1, 0, 0, 1);
}

vec4 fetch(sampler2D tex, int index)
{
	int width = textureSize(tex, 0).x;
	int y = index / width;
	int x = index % width;
	return texelFetch(tex, ivec2(x, y), 0);
}

float intersectTriangle(vec3 v1, vec3 v2, vec3 v3, vec3 origin, vec3 direction)
{
	vec3 vertices[3];
	vertices[0] = v1;
	vertices[1] = v2;
	vertices[2] = v3;
	vec3 lines[3], r2v, n;
	int mask_pos = 0;
	int mask_neg = 0;
	for(int i=0; i<3; i++)
	{
		lines[i] = vertices[(i+1)%3] - vertices[i];
		r2v = vertices[i] - origin;
		n = cross(lines[i], r2v);
		mask_pos += (dot(direction, n) >= 0) << i;
		mask_neg += (dot(direction, n) <= 0) << i;
	}
	n = cross(lines[0], lines[1]);
	n = normalize(n);
	float orthoDist = dot(r2v, n);
	float d = dot(direction, n);
	if(d == 0 || orthoDist/d < 0 || (mask_pos!=7 && mask_neg != 7))
		return -1;
	float dist = orthoDist/d;
	return dist;
}
