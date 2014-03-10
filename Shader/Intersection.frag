#version 150
#extension GL_EXT_gpu_shader4 : enable

uniform int nQueries;

uniform sampler2D queryRayOriginTex;
uniform sampler2D queryRayDirectionTex;

uniform sampler2D kdTreeNodes_pi_ri_ls_le;
uniform sampler2D kdTreeNodes_minCoords;
uniform sampler2D kdTreeNodes_maxCoords;
uniform sampler2D leafVertexPosition1;
uniform sampler2D leafVertexPosition2;
uniform sampler2D leafVertexPosition3;



float intersectTriangle(vec3 v1, vec3 v2, vec3 v3, vec3 origin, vec3 direction);
float intersectBBox(vec3 minCoord, vec3 maxCoord, vec3 origin, vec3 direction);

vec4 fetch(sampler2D tex, int index);

vec2 min_dist_tid(vec2 dt1, vec2 dt2);

void main()
{
	//gl_FragColor = vec4(textureSize(kdTreeNodes_pi_ri_ls_le, 0).x, textureSize(kdTreeNodes_pi_ri_ls_le, 0).y, 0, 0);
	//return;
	
	int queryIndex = int(gl_FragCoord.y) * textureSize(queryRayOriginTex, 0).x + int(gl_FragCoord.x);
	if(queryIndex >= nQueries)
	{
		gl_FragColor = vec4(0);
		return;
	}
	vec4 origin = fetch(queryRayOriginTex, queryIndex);
	vec3 direction = fetch(queryRayDirectionTex, queryIndex).xyz;
	
	vec2 dist_tid = vec2(-1, -1);
	
	int stackSize = 0;
	int stack[50];
	int nodeIndex = 0;
	
	while(true)
	{
		float bbDist = intersectBBox(fetch(kdTreeNodes_minCoords, nodeIndex).xyz,
			fetch(kdTreeNodes_maxCoords, nodeIndex).xyz, origin.xyz, direction);
		
		if(bbDist < 0)
		{
			if(stackSize == 0)
				break;
			stackSize--;
			nodeIndex = stack[stackSize]; // pop
		}
		vec4 node = fetch(kdTreeNodes_pi_ri_ls_le, nodeIndex);
		if(node[1] < 0) // is leaf
		{
			for(int ti=int(node[2]); ti<int(node[3]); ti++)
			{
				vec4 v1 = fetch(leafVertexPosition1, ti);
				if(int(origin.w) == int(v1.w))
					continue;
				vec3 v2 = fetch(leafVertexPosition2, ti).xyz;
				vec3 v3 = fetch(leafVertexPosition3, ti).xyz;
				float dist = intersectTriangle(v1.xyz, v2, v3, origin.xyz, direction);
				dist_tid = min_dist_tid(dist_tid, vec2(dist, v1.w));
			}
			if(stackSize == 0)
				break;
			stackSize--;
			nodeIndex = stack[stackSize]; // pop
		}
		else
		{
			stack[stackSize] = int(node[1]);
			stackSize++; // push
			nodeIndex++; // go left
		}
	}
	
	gl_FragColor = vec4(dist_tid, 0, 0);
}

vec2 min_dist_tid(vec2 dt1, vec2 dt2)
{
	if(dt1.x < 0)
		return dt2;
	if(dt2.x < 0)
		return dt1;
	return dt1.x > dt2.x ? dt2 : dt1;
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
		mask_pos += (int(dot(direction, n) >= 0)) << i;
		mask_neg += (int(dot(direction, n) <= 0)) << i;
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

float intersectBBox(vec3 minCoord, vec3 maxCoord, vec3 origin, vec3 direction)
{
	vec3 diff = maxCoord - minCoord;
	vec3 intPoint;
	float minDist = -1;
	float coord[6];
	coord[0] = minCoord[0];
	coord[1] = minCoord[1];
	coord[2] = minCoord[2];
	coord[3] = maxCoord[0];
	coord[4] = maxCoord[1];
	coord[5] = maxCoord[2];
	
	minCoord -= diff*(1e-5);
	maxCoord += diff*(1e-5);
	for(int i=0; i<6; i++)
	{
		if(direction[i%3] == 0)
			continue;
		
		float dist = (coord[i] - origin[i%3]) / direction[i%3];
		intPoint = origin + dist*direction;
		bool inter = dist >= 0;
		inter = inter && minCoord[(i+1)%3] <= intPoint[(i+1)%3];
		inter = inter && maxCoord[(i+1)%3] >= intPoint[(i+1)%3];
		inter = inter && minCoord[(i+2)%3] <= intPoint[(i+2)%3];
		inter = inter && maxCoord[(i+2)%3] >= intPoint[(i+2)%3];
		minDist = (inter && (dist < minDist || minDist < 0)) ? dist : minDist;
	}
	return minDist;
}
