using System;
using System.Collections.Generic;
using System.Linq;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using UnityEngine;

public class MeshBuilderDualContouringI : Singleton<MeshBuilderDualContouringI>
{
    [Tooltip("Value from which the vertices are inside the figure")]
    [Range(0, 255)]
    public int isoLevel = 128;

    [Tooltip("Allow to get a middle point between the voxel vertices in function of the weight of the vertices")]
    public bool interpolate = false;

    public bool dualContouring = true;

    public bool useDeprecated;

    public float terrainSurface = 0.5f;

    /// <summary>
    /// Method that calculate cubes, vertex and mesh in that order of a chunk.
    /// </summary>
    /// <param name="b"> data of the chunk</param>
    public Mesh BuildChunk(byte[] b)
    {
        if (useDeprecated)
            return BuildChunkDeprecated(b);

        var buildChunkJob = new BuildChunkJob
        {
            chunkData = new NativeArray<byte>(b, Allocator.TempJob),
            isoLevel = isoLevel,
            interpolate = interpolate,
            vertex = new NativeList<float3>(500, Allocator.TempJob),
            triangles = new NativeList<int>(500, Allocator.TempJob),
            uv = new NativeList<float2>(100, Allocator.TempJob),
            dualContouring = dualContouring,
            terrainSurface = 0.5f
        };
        var jobHandle = buildChunkJob.Schedule();
        jobHandle.Complete();

        //Get all the data from the jobs and use to generate a Mesh
        var meshGenerated = new Mesh();

        var meshVert = new Vector3[buildChunkJob.vertex.Length];

        for (var i = 0; i < buildChunkJob.vertex.Length; i++)
            meshVert[i] = buildChunkJob.vertex[i];
        meshGenerated.vertices = meshVert;

        var meshTriangles = new int[buildChunkJob.triangles.Length];
        for (var i = 0; i < buildChunkJob.triangles.Length; i++)
            meshTriangles[i] = i;
        meshGenerated.triangles = meshTriangles;

        var meshUV = new Vector2[buildChunkJob.vertex.Length];
        for (var i = 0; i < buildChunkJob.vertex.Length; i++)
            meshUV[i] = buildChunkJob.uv[i];
        meshGenerated.uv = meshUV;

        meshGenerated.RecalculateNormals();
        meshGenerated.RecalculateTangents();

        //Dispose (Clear the jobs NativeLists)
        buildChunkJob.vertex.Dispose();
        buildChunkJob.triangles.Dispose();
        buildChunkJob.uv.Dispose();
        buildChunkJob.chunkData.Dispose();

        return meshGenerated;
    }

    //This old code was adapted in the "BuildChunkJob" script and don't used anymore. (Stay if someone want to use the )

    #region Original code (Deprecated)

    /// <summary>
    /// Method that calculate cubes, vertex and mesh in that order of a chunk.
    /// </summary>
    /// <param name="b"> data of the chunk</param>
    public Mesh BuildChunkDeprecated(byte[] b)
    {
        var vertex = new List<Vector3>();
        var triangles = new List<int>();
        var matVert = new List<Vector2>();

        for (var y = 0; y < Constants.MAX_HEIGHT; y++)//height
        {
            for (var z = 1; z < Constants.CHUNK_SIZE + 1; z++)//column, start at 1, because Z axis is inverted and need -1 as offset
            {
                for (var x = 0; x < Constants.CHUNK_SIZE; x++)//line
                {
                    var cube = new Vector4[8];
                    var mat = Constants.NUMBER_MATERIALS;
                    cube[0] = CalculateVertexChunk(x, y, z, b, ref mat);
                    cube[1] = CalculateVertexChunk(x + 1, y, z, b, ref mat);
                    cube[2] = CalculateVertexChunk(x + 1, y, z - 1, b, ref mat);
                    cube[3] = CalculateVertexChunk(x, y, z - 1, b, ref mat);
                    cube[4] = CalculateVertexChunk(x, y + 1, z, b, ref mat);
                    cube[5] = CalculateVertexChunk(x + 1, y + 1, z, b, ref mat);
                    cube[6] = CalculateVertexChunk(x + 1, y + 1, z - 1, b, ref mat);
                    cube[7] = CalculateVertexChunk(x, y + 1, z - 1, b, ref mat);

                    CalculateVertex(cube, mat, matVert, vertex, triangles);
                }
            }
        }
        return BuildMesh(vertex, triangles, matVert);
    }

    /// <summary>
    /// It generate a mesh from a group of vertex. Flat shading type.(Deprecated)
    /// </summary>
    public Mesh BuildMesh(List<Vector3> vertex, List<int> triangles, List<Vector2> textures = null)
    {
        var mesh = new Mesh();

        mesh.vertices = vertex.ToArray();
        mesh.triangles = triangles.ToArray();

        if (textures != null)
            mesh.uv = textures.ToArray();

        return mesh;
    }

    /// <summary>
    ///  Calculate the vertices of the voxels, get the vertices of the triangulation table and his position in the world. Also check materials of that vertex (UV position).(Deprecated)
    /// </summary>
    public void CalculateVertex(Vector4[] cube, int colorVert, List<Vector2> matVert, List<Vector3> vertex, List<int> triangles)
    {
        //Values above isoLevel are inside the figure, value of 0 means that the cube is entirely inside of the figure.
        var cubeindex = 0;
        if (cube[0].w < isoLevel) cubeindex |= 1;
        if (cube[1].w < isoLevel) cubeindex |= 2;
        if (cube[2].w < isoLevel) cubeindex |= 4;
        if (cube[3].w < isoLevel) cubeindex |= 8;
        if (cube[4].w < isoLevel) cubeindex |= 16;
        if (cube[5].w < isoLevel) cubeindex |= 32;
        if (cube[6].w < isoLevel) cubeindex |= 64;
        if (cube[7].w < isoLevel) cubeindex |= 128;

        //List<Vector3> vertexArray = new List<Vector3>();

        for (var i = 0; Constants.triTable[cubeindex, i] != -1; i++)
        {
            var v1 = Constants.cornerIndexAFromEdge[Constants.triTable[cubeindex, i]];
            var v2 = Constants.cornerIndexBFromEdge[Constants.triTable[cubeindex, i]];

            var indice = vertex.Count;

            if (dualContouring)
            {
                try
                {
                    // Get the terrain values at either end of our current edge from the cube array created above.
                    var vert1Sample = cube[EdgeIndexes[(indice / 2) % 24]].w;
                    var vert2Sample = cube[EdgeIndexes[(indice / 2 + 1) % 24]].w;

                    // Calculate the difference between the terrain values.
                    var difference = vert2Sample - vert1Sample;

                    // If the difference is 0, then the terrain passes through the middle.
                    if (difference == 0)
                        difference = terrainSurface;
                    else
                        difference = (terrainSurface - vert1Sample) / difference;

                    // Calculate the point along the edge that passes through.

                    var p1 = cube[v1];
                    var p2 = cube[v2];

                    var vertPosition = p1 + ((p2 - p1) * difference);
                    triangles.Add(VertForIndice(vertPosition, vertex));
                }
                catch (Exception ex)
                {
                    //Debug.Log(indice / 2 + 1);
                    Debug.LogException(ex);
                }
            }
            else
            {
                //vertexArray.Add(MiddlePointVertex(cube[v1], cube[v2]));

                if (interpolate)
                    vertex.Add(InterpolateVertex(cube[v1], cube[v2], cube[v1].w, cube[v2].w));
                else
                    vertex.Add(MiddlePointVertex(cube[v1], cube[v2]));

                triangles.Add(indice);
            }

            const float uvOffset = 0.01f; //Small offset for avoid pick pixels of other textures
            //NEED REWORKING FOR CORRECT WORKING, now have problems with the directions of the uv
            if (i % 6 == 0)
                matVert.Add(new Vector2(Constants.MATERIAL_SIZE * (colorVert % Constants.MATERIAL_FOR_ROW) + Constants.MATERIAL_SIZE - uvOffset,
                                  1 - Constants.MATERIAL_SIZE * Mathf.Floor(colorVert / Constants.MATERIAL_FOR_ROW) - uvOffset));
            else if (i % 6 == 1)
                matVert.Add(new Vector2(Constants.MATERIAL_SIZE * (colorVert % Constants.MATERIAL_FOR_ROW) + Constants.MATERIAL_SIZE - uvOffset,
                                  1 - Constants.MATERIAL_SIZE * Mathf.Floor(colorVert / Constants.MATERIAL_FOR_ROW) - Constants.MATERIAL_SIZE + uvOffset));
            else if (i % 6 == 2)
                matVert.Add(new Vector2(Constants.MATERIAL_SIZE * (colorVert % Constants.MATERIAL_FOR_ROW) + uvOffset,
                                  1 - Constants.MATERIAL_SIZE * Mathf.Floor(colorVert / Constants.MATERIAL_FOR_ROW) - uvOffset));
            else if (i % 6 == 3)
                matVert.Add(new Vector2(Constants.MATERIAL_SIZE * (colorVert % Constants.MATERIAL_FOR_ROW) + Constants.MATERIAL_SIZE - uvOffset,
                                  1 - Constants.MATERIAL_SIZE * Mathf.Floor(colorVert / Constants.MATERIAL_FOR_ROW) - Constants.MATERIAL_SIZE + uvOffset));
            else if (i % 6 == 4)
                matVert.Add(new Vector2(Constants.MATERIAL_SIZE * (colorVert % Constants.MATERIAL_FOR_ROW) + uvOffset,
                                   1 - Constants.MATERIAL_SIZE * Mathf.Floor(colorVert / Constants.MATERIAL_FOR_ROW) - Constants.MATERIAL_SIZE + uvOffset));
            else if (i % 6 == 5)
                matVert.Add(new Vector2(Constants.MATERIAL_SIZE * (colorVert % Constants.MATERIAL_FOR_ROW + uvOffset),
                                  1 - Constants.MATERIAL_SIZE * Mathf.Floor(colorVert / Constants.MATERIAL_FOR_ROW) - uvOffset));
        }
    }

    public int VertForIndice(Vector3 vert, List<Vector3> vertex)
    {
        // Loop through all the vertices currently in the vertices list.
        for (var i = 0; i < vertex.Count; i++)
        {
            // If we find a vert that matches ours, then simply return this index.
            if (vertex[i].Equals(vert))
                return i;
        }

        // If we didn't find a match, add this vert to the list and return last index.
        vertex.Add(vert);
        return vertex.Count - 1;
    }

    /// <summary>
    /// Calculate the data of a vertex of a voxel.(Deprecated)
    /// </summary>
    private Vector4 CalculateVertexChunk(int x, int y, int z, byte[] b, ref int colorVoxel)
    {
        var index = (x + z * Constants.CHUNK_VERTEX_SIZE + y * Constants.CHUNK_VERTEX_AREA) * Constants.CHUNK_POINT_BYTE;
        int material = b[index + 1];
        if (b[index] >= isoLevel && material < colorVoxel)
            colorVoxel = material;
        return new Vector4(
            (x - Constants.CHUNK_SIZE / 2) * Constants.VOXEL_SIDE,
            (y - Constants.MAX_HEIGHT / 2) * Constants.VOXEL_SIDE,
            (z - Constants.CHUNK_SIZE / 2) * Constants.VOXEL_SIDE,
            b[index]);
    }

    /// <summary>
    /// Overload of the CalculateVertex method but without material calculations.
    /// </summary>
    public List<Vector3> CalculateVertex(Vector4[] cube)
    {
        //Values above isoLevel are inside the figure, value of 0 means that the cube is entirely inside of the figure.(Deprecated)
        var cubeindex = 0;
        if (cube[0].w < isoLevel) cubeindex |= 1;
        if (cube[1].w < isoLevel) cubeindex |= 2;
        if (cube[2].w < isoLevel) cubeindex |= 4;
        if (cube[3].w < isoLevel) cubeindex |= 8;
        if (cube[4].w < isoLevel) cubeindex |= 16;
        if (cube[5].w < isoLevel) cubeindex |= 32;
        if (cube[6].w < isoLevel) cubeindex |= 64;
        if (cube[7].w < isoLevel) cubeindex |= 128;

        var vertexArray = new List<Vector3>();

        for (var i = 0; Constants.triTable[cubeindex, i] != -1; i++)
        {
            var v1 = Constants.cornerIndexAFromEdge[Constants.triTable[cubeindex, i]];
            var v2 = Constants.cornerIndexBFromEdge[Constants.triTable[cubeindex, i]];

            if (interpolate)
                vertexArray.Add(InterpolateVertex(cube[v1], cube[v2], cube[v1].w, cube[v2].w));
            else
                vertexArray.Add(MiddlePointVertex(cube[v1], cube[v2]));
        }

        return vertexArray;
    }

    //HelpMethods

    /// <summary>
    /// Calculate a point between two vertex using the weight of each vertex , used in interpolation voxel building.(Deprecated)
    /// </summary>
    public Vector3 InterpolateVertex(Vector3 p1, Vector3 p2, float val1, float val2)
    {
        return Vector3.Lerp(p1, p2, (isoLevel - val1) / (val2 - val1));
    }

    /// <summary>
    /// Calculate the middle point between two vertex, for no interpolation voxel building.(Deprecated)
    /// </summary>
    public Vector3 MiddlePointVertex(Vector3 p1, Vector3 p2)
    {
        return (p1 + p2) / 2;
    }

    public static readonly int[] EdgeIndexes = new int[] {
        0, 1,
        1, 2,
        3, 2,
        0, 3,
        4, 5,
        5, 6,
        7, 6,
        4, 7,
        0, 4,
        1, 5,
        2, 6,
        3, 7
    };

    #endregion Original code (Deprecated)
}