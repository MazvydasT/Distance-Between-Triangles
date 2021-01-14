using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Windows.Media.Media3D;

using Geometry3D;
using Geometry3DTests;

namespace Geometry3DTestsDebugger
{
    class Program
    {
        static void Main(string[] args)
        {
            //new GeometryTests().Point2Point();
            //new GeometryTests().Point2Line();
            //new GeometryTests().Point2Triangle();
            
            //new GeometryTests().Triangle2Triangle();

            //GetMinBoundingSphere(new Triangle3D(new Point3D(-136.031723022, 101.843025208, 0), new Point3D(-136.061553955, 16.804965973, 0), new Point3D(-72.061553955, 16.782518387, 0)));
            //GetMinBoundingSphere(new Triangle3D(new Point3D(0.01, 0.01, 0.01), new Point3D(2.331, 0.01, 0.01), new Point3D(0.01, 112.459871, 0.01)));
            //var asd = GetMinBoundingSphere(new Triangle3D(new Point3D(113.92, 40.024, 0), new Point3D(102.424, 8.438, 0), new Point3D(140.424, 8.438, 0)));

            //return;

            var data = GeometryTestData.Triangle2TriangleData;
            var dataLength = data.Length;

            var uniquePointTracker1 = new Dictionary<Point3D, int>(dataLength * 3);
            var uniquePointTracker2 = new Dictionary<Point3D, int>(dataLength * 3);

            var vertexBuffer1 = new List<Point3D>();
            var indexBuffer1 = new List<int>();

            var vertexBuffer2 = new List<Point3D>();
            var indexBuffer2 = new List<int>();

            //var mesh1Data = new Triangle3D[dataLength];
            //var mesh2Data = new Triangle3D[dataLength];

            for (int i = 0; i < dataLength; ++i)
            {
                var item = data[i];
                var g1 = (Triangle3D)item.g1;
                var g2 = (Triangle3D)item.g2;

                int vertexIndex;

                if (!uniquePointTracker1.TryGetValue(g1.point1, out vertexIndex))
                {
                    vertexIndex = indexBuffer1.Count;
                    uniquePointTracker1.Add(g1.point1, vertexIndex);

                    vertexBuffer1.Add(g1.point1);
                }

                indexBuffer1.Add(vertexIndex);

                if (!uniquePointTracker1.TryGetValue(g1.point2, out vertexIndex))
                {
                    vertexIndex = indexBuffer1.Count;
                    uniquePointTracker1.Add(g1.point2, vertexIndex);

                    vertexBuffer1.Add(g1.point2);
                }

                indexBuffer1.Add(vertexIndex);

                if (!uniquePointTracker1.TryGetValue(g1.point3, out vertexIndex))
                {
                    vertexIndex = indexBuffer1.Count;
                    uniquePointTracker1.Add(g1.point3, vertexIndex);

                    vertexBuffer1.Add(g1.point3);
                }

                indexBuffer1.Add(vertexIndex);

                //############################################################

                if (!uniquePointTracker2.TryGetValue(g2.point1, out vertexIndex))
                {
                    vertexIndex = indexBuffer2.Count;
                    uniquePointTracker2.Add(g2.point1, vertexIndex);

                    vertexBuffer2.Add(g2.point1);
                }

                indexBuffer2.Add(vertexIndex);

                if (!uniquePointTracker2.TryGetValue(g2.point2, out vertexIndex))
                {
                    vertexIndex = indexBuffer2.Count;
                    uniquePointTracker2.Add(g2.point2, vertexIndex);

                    vertexBuffer2.Add(g2.point2);
                }

                indexBuffer2.Add(vertexIndex);

                if (!uniquePointTracker2.TryGetValue(g2.point3, out vertexIndex))
                {
                    vertexIndex = indexBuffer2.Count;
                    uniquePointTracker2.Add(g2.point3, vertexIndex);

                    vertexBuffer2.Add(g2.point3);
                }

                indexBuffer2.Add(vertexIndex);
            }

            double shortestDistanceBetweenMeshes = -1;
            LineSegment3D shortestLineBetweenMeshes = null;

            var sw = new Stopwatch();
            sw.Start();

            var t2s = new Triangle3D[indexBuffer2.Count / 3];
            var t2BoundingCircles = new Sphere3D[t2s.Length];

            for (int i1 = 0, c1 = indexBuffer1.Count; i1 < c1; i1 += 3)
            {
                var t1 = new Triangle3D(vertexBuffer1[indexBuffer1[i1]], vertexBuffer1[indexBuffer1[i1 + 1]], vertexBuffer1[indexBuffer1[i1 + 2]]);

                //var boundingCircle1 = GetBoundingSphere(t1);
                var boundingCircle1 = Utilities.GetMinimumBoundingSphere(t1);
                var boundingCircleCentre = boundingCircle1.centrePoint;
                var boundingCircleRadius = boundingCircle1.radius;

                for (int i2 = 0, c2 = indexBuffer2.Count; i2 < c2; i2 += 3)
                {
                    Triangle3D t2;
                    Sphere3D boundingCircle2;

                    var t2Index = i2/3;

                    if (i1 == 0)
                    {
                        t2 = new Triangle3D(vertexBuffer2[indexBuffer2[i2]], vertexBuffer2[indexBuffer2[i2 + 1]], vertexBuffer2[indexBuffer2[i2 + 2]]);
                        //boundingCircle2 = GetBoundingSphere(t2);
                        boundingCircle2 = Utilities.GetMinimumBoundingSphere(t2);

                        t2s[t2Index] = t2;
                        t2BoundingCircles[t2Index] = boundingCircle2;
                    }

                    else
                    {
                        t2 = t2s[t2Index];
                        boundingCircle2 = t2BoundingCircles[t2Index];
                    }

                    var distanceBetweenBoundingSpheres = Math.Max((boundingCircleCentre - boundingCircle2.centrePoint).Length - boundingCircleRadius - boundingCircle2.radius, 0);

                    if (shortestDistanceBetweenMeshes == -1 || distanceBetweenBoundingSpheres < shortestDistanceBetweenMeshes)
                    {
                        var shortestLineBetweenTriangles = Utilities.ShortestLineSegmentBetween(t1, t2);
                        var shortestDistanceBetweenTriangles = shortestLineBetweenTriangles.Vector.Length;

                        if (shortestDistanceBetweenMeshes == -1 || shortestDistanceBetweenTriangles < shortestDistanceBetweenMeshes)
                        {
                            shortestDistanceBetweenMeshes = shortestDistanceBetweenTriangles;
                            shortestLineBetweenMeshes = shortestLineBetweenTriangles;
                        }
                    }

                    /*var shortestLineBetweenTriangles = Utilities.ShortestLineSegmentBetween(t1, t2);
                    var shortestDistanceBetweenTriangles = shortestLineBetweenTriangles.Vector.Length;

                    if (shortestDistanceBetweenMeshes == -1 || shortestDistanceBetweenTriangles < shortestDistanceBetweenMeshes)
                    {
                        shortestDistanceBetweenMeshes = shortestDistanceBetweenTriangles;
                        shortestLineBetweenMeshes = shortestLineBetweenTriangles;
                    }*/
                }
            }

            /*for(int i1 = 0, c1 = mesh1Data.Length; i1 < c1; ++i1)
            {
                var t1 = mesh1Data[i1];

                for (int i2 = 0, c2 = mesh2Data.Length; i2 < c2; ++i2)
                {
                    var t2 = mesh2Data[i2];

                    var shortestLineBetweenTriangles = Utilities.ShortestLineSegmentBetween(t1, t2);
                    var shortestDistanceBetweenTriangles = shortestLineBetweenTriangles.Vector.Length;

                    if (i1 == 0 || shortestDistanceBetweenTriangles < shortestDistanceBetweenMeshes)
                    {
                        shortestDistanceBetweenMeshes = shortestDistanceBetweenTriangles;
                        shortestLineBetweenMeshes = shortestLineBetweenTriangles;
                    }
                }
            }*/

            sw.Stop();

            var el = sw.ElapsedMilliseconds;
            Console.WriteLine(el);
        }
    }
}
