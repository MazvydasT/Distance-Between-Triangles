using System;
using System.Windows.Media.Media3D;

namespace Geometry3D
{
    public class LineSegment3D
    {
        public readonly Point3D point1;
        public readonly Point3D point2;

        public Vector3D Vector { get { return point2 - point1; } }

        public LineSegment3D(Point3D point1, Point3D point2)
        {
            this.point1 = point1;
            this.point2 = point2;
        }
    }

    public class Triangle3D
    {
        public readonly Point3D point1;
        public readonly Point3D point2;
        public readonly Point3D point3;

        public LineSegment3D Edge1 { get { return new LineSegment3D(point1, point2); } }
        public LineSegment3D Edge2 { get { return new LineSegment3D(point2, point3); } }
        public LineSegment3D Edge3 { get { return new LineSegment3D(point3, point1); } }

        public Vector3D Normal { get { return Vector3D.CrossProduct(Edge1.Vector, Edge2.Vector); } }

        public Triangle3D(Point3D point1, Point3D point2, Point3D point3)
        {
            this.point1 = point1;
            this.point2 = point2;
            this.point3 = point3;
        }
    }

    public class Sphere3D
    {
        public readonly Point3D centrePoint;
        public readonly double radius;

        public Sphere3D(Point3D centrePoint, double radius)
        {
            this.centrePoint = centrePoint;
            this.radius = radius;
        }
    }

    public static class Utilities
    {
        public static LineSegment3D GetShortest(Func<LineSegment3D>[] getLineSegmentFunctions)
        {
            LineSegment3D shortestDistancePoints = null;
            double shortestDistanceSquared = 0;

            for (int i = 0, c = getLineSegmentFunctions.Length; i < c; ++i)
            {
                var lineSegment = getLineSegmentFunctions[i]();
                var lineVector = lineSegment.Vector;
                var lengthSquared = lineVector.LengthSquared;

                if (lengthSquared == 0)
                    return lineSegment;

                if (i == 0 || lengthSquared < shortestDistanceSquared)
                {
                    shortestDistancePoints = lineSegment;
                    shortestDistanceSquared = lengthSquared;
                }
            }

            return shortestDistancePoints;
        }

        public static Sphere3D GetMinimumBoundingSphere(Triangle3D triangle)
        {
            var point1 = triangle.point1;
            var point2 = triangle.point2;
            var point3 = triangle.point3;

            var vector12 = point2 - point1;
            var vector23 = point3 - point2;
            var vector31 = point1 - point3;

            var line12LengthSquared = vector12.LengthSquared;
            var line23LengthSquared = vector23.LengthSquared;
            var line31LengthSquared = vector31.LengthSquared;

            Point3D longestEdgePoint1, longestEdgePoint2;
            double longestEdgeLengthSquared, otherEdge1LengthSquared, otherEdge2LengthSquared;

            if (line12LengthSquared > line23LengthSquared && line12LengthSquared > line31LengthSquared)
            {
                longestEdgePoint1 = point1;
                longestEdgePoint2 = point2;
                longestEdgeLengthSquared = line12LengthSquared;
                otherEdge1LengthSquared = line23LengthSquared;
                otherEdge2LengthSquared = line31LengthSquared;
            }

            else if (line23LengthSquared > line12LengthSquared && line23LengthSquared > line31LengthSquared)
            {
                longestEdgePoint1 = point2;
                longestEdgePoint2 = point3;
                longestEdgeLengthSquared = line23LengthSquared;
                otherEdge1LengthSquared = line12LengthSquared;
                otherEdge2LengthSquared = line31LengthSquared;
            }

            else
            {
                longestEdgePoint1 = point3;
                longestEdgePoint2 = point1;
                longestEdgeLengthSquared = line31LengthSquared;
                otherEdge1LengthSquared = line12LengthSquared;
                otherEdge2LengthSquared = line23LengthSquared;
            }

            // workout longest edge

            Point3D centrePoint;
            double radius;

            if (longestEdgeLengthSquared < otherEdge1LengthSquared + otherEdge2LengthSquared)
            {
                // Triangle is acute

                var triangleNormal = Vector3D.CrossProduct(vector12, vector23);

                var bisector1Vector = Vector3D.CrossProduct(triangleNormal, vector12);
                var bisector2Vector = Vector3D.CrossProduct(triangleNormal, vector23);

                var midPoint1 = new Point3D((point1.X + point2.X) / 2, (point1.Y + point2.Y) / 2, (point1.Z + point2.Z) / 2);
                var midPoint2 = new Point3D((point2.X + point3.X) / 2, (point2.Y + point3.Y) / 2, (point2.Z + point3.Z) / 2);

                var normalBetweenBisectorAndTriangleNormal = Vector3D.CrossProduct(bisector1Vector, triangleNormal);

                centrePoint = midPoint2 + (Vector3D.DotProduct(midPoint1 - midPoint2, normalBetweenBisectorAndTriangleNormal) / Vector3D.DotProduct(bisector2Vector, normalBetweenBisectorAndTriangleNormal)) * bisector2Vector;
                radius = (centrePoint - point1).Length;
            }

            else
            {
                // Triangle is obtuse or right

                var midPointX = (longestEdgePoint1.X + longestEdgePoint2.X) / 2;
                var midPointY = (longestEdgePoint1.Y + longestEdgePoint2.Y) / 2;
                var midPointZ = (longestEdgePoint1.Z + longestEdgePoint2.Z) / 2;

                centrePoint = new Point3D(midPointX, midPointY, midPointZ);
                radius = Math.Sqrt(longestEdgeLengthSquared) / 2;
            }

            return new Sphere3D(centrePoint, radius);
        }

        public static Point3D Project(Point3D point, LineSegment3D lineSegment)
        {
            var point1ToPointVector = point - lineSegment.point1;
            var point1ToPoint2Vector = lineSegment.point2 - lineSegment.point1;

            return lineSegment.point1 + Vector3D.DotProduct(point1ToPointVector, point1ToPoint2Vector) / Vector3D.DotProduct(point1ToPoint2Vector, point1ToPoint2Vector) * point1ToPoint2Vector;
        }

        public static Point3D Project(Point3D point, Triangle3D triangle)
        {
            var linePerpendicularToTriangle = new LineSegment3D(triangle.point1 + triangle.Normal, triangle.point1);
            var pointProjectedToLinePerpendicularToTriangle = Utilities.Project(point, linePerpendicularToTriangle);
            var point1ToProjectedPointVector = pointProjectedToLinePerpendicularToTriangle - triangle.point1;

            return point - point1ToProjectedPointVector;
        }

        public static bool AreIncident(Point3D point, LineSegment3D lineSegment)
        {
            return Math.Abs(((lineSegment.point1 - point).Length + (lineSegment.point2 - point).Length) - lineSegment.Vector.Length) < 0.001;
        }
        
        public static bool AreIncident(Point3D point, Triangle3D triangle)
        {
            var v0 = triangle.point3 - triangle.point1;
            var v1 = triangle.point2 - triangle.point1;
            var v2 = point - triangle.point1;

            var dot00 = Vector3D.DotProduct(v0, v0);
            var dot01 = Vector3D.DotProduct(v0, v1);
            var dot02 = Vector3D.DotProduct(v0, v2);
            var dot11 = Vector3D.DotProduct(v1, v1);
            var dot12 = Vector3D.DotProduct(v1, v2);

            var invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
            var u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            var v = (dot00 * dot12 - dot01 * dot02) * invDenom;

            return (u >= 0) && (v >= 0) && (u + v < 1);
        }
        
        public static bool AreParallel(LineSegment3D lineSegment1, LineSegment3D lineSegment2)
        {
            var point1 = lineSegment1.point1;
            var point2 = lineSegment1.point2;

            var point1ProjectedTolineSegment2 = Utilities.Project(point1, lineSegment2);
            var point2ProjectedToLineSegment2 = Utilities.Project(point2, lineSegment2);

            return Utilities.ShortestLineSegmentBetween(point2ProjectedToLineSegment2, point2ProjectedToLineSegment2).Vector.Length == lineSegment1.Vector.Length;
        }
        
        public static bool AreParallel(LineSegment3D lineSegment, Triangle3D triangle)
        {
            return Vector3D.DotProduct(lineSegment.Vector, triangle.Normal) == 0;
        }
        
        public static LineSegment3D ShortestLineSegmentBetween(Point3D point1, Point3D point2)
        {
            return new LineSegment3D(point1, point2);
        }
        
        public static LineSegment3D ShortestLineSegmentBetween(Point3D point, LineSegment3D lineSegment)
        {
            var projectedPoint = Utilities.Project(point, lineSegment);

            if (Utilities.AreIncident(projectedPoint, lineSegment))
                return Utilities.ShortestLineSegmentBetween(point, projectedPoint);

            return Utilities.GetShortest(new Func<LineSegment3D>[]
            {
                () => Utilities.ShortestLineSegmentBetween(point, lineSegment.point1),
                () => Utilities.ShortestLineSegmentBetween(point, lineSegment.point2)
            });
        }

        public static LineSegment3D ShortestLineSegmentBetween(Point3D point, Triangle3D triangle)
        {
            var pointProjectedToTrianglePlane = Utilities.Project(point, triangle);

            if (Utilities.AreIncident(pointProjectedToTrianglePlane, triangle))
                return Utilities.ShortestLineSegmentBetween(point, pointProjectedToTrianglePlane);
            
            else
                return Utilities.GetShortest(new Func<LineSegment3D>[]
                {
                    () => Utilities.ShortestLineSegmentBetween(point, triangle.Edge1),
                    () => Utilities.ShortestLineSegmentBetween(point, triangle.Edge2),
                    () => Utilities.ShortestLineSegmentBetween(point, triangle.Edge3)
                });
        }
        
        public static LineSegment3D ShortestLineSegmentBetween(LineSegment3D lineSegment1, LineSegment3D lineSegment2)
        {
            var lineSegment1Point1 = lineSegment1.point1;
            var lineSegment1Point2 = lineSegment1.point2;

            var lineSegment2Point1 = lineSegment2.point1;
            var lineSegment2Point2 = lineSegment2.point2;

            if (Utilities.AreParallel(lineSegment1, lineSegment2))
                return Utilities.GetShortest(new Func<LineSegment3D>[]
                {
                    () => Utilities.ShortestLineSegmentBetween(lineSegment1Point1, lineSegment2),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment1Point2, lineSegment2),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment2Point1, lineSegment1),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment2Point2, lineSegment1)
                });

            var normalVector = Vector3D.CrossProduct(lineSegment1.Vector, lineSegment2.Vector);
            
            var lineSegment1Vector = lineSegment1.Vector;
            var lineSegment2Vector = lineSegment2.Vector;

            var n1 = Vector3D.CrossProduct(lineSegment1Vector, normalVector);
            var n2 = Vector3D.CrossProduct(lineSegment2Vector, normalVector);

            var pointOnLine1 = lineSegment1Point1 + (Vector3D.DotProduct(lineSegment2Point1 - lineSegment1Point1, n2) / Vector3D.DotProduct(lineSegment1Vector, n2)) * lineSegment1Vector;
            var pointOnLine2 = lineSegment2Point1 + (Vector3D.DotProduct(lineSegment1Point1 - lineSegment2Point1, n1) / Vector3D.DotProduct(lineSegment2Vector, n1)) * lineSegment2Vector;

            if (!Utilities.AreIncident(pointOnLine1, lineSegment1) || !Utilities.AreIncident(pointOnLine2, lineSegment2))
                return Utilities.GetShortest(new Func<LineSegment3D>[]
                {
                    () => Utilities.ShortestLineSegmentBetween(lineSegment1Point1, lineSegment2),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment1Point2, lineSegment2),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment2Point1, lineSegment1),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment2Point2, lineSegment1)
                });

            return Utilities.ShortestLineSegmentBetween(pointOnLine1, pointOnLine2);
        }
        
        public static LineSegment3D ShortestLineSegmentBetween(LineSegment3D lineSegment, Triangle3D triangle)
        {
            var triangleEdge1 = triangle.Edge1;
            var triangleEdge2 = triangle.Edge2;
            var triangleEdge3 = triangle.Edge3;

            if (Utilities.AreParallel(lineSegment, triangle))
            {
                return Utilities.GetShortest(new Func<LineSegment3D>[]
                {
                    () => Utilities.ShortestLineSegmentBetween(lineSegment.point1, triangle),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment.point2, triangle),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge1),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge2),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge3)
                });
            }

            var triangleNormal = triangle.Normal;
            var lineSegmentVector = lineSegment.Vector;

            var linePlaneIntersectionPoint = (Vector3D.DotProduct(triangle.point1 - lineSegment.point1, triangleNormal) / Vector3D.DotProduct(lineSegmentVector, triangleNormal)) * lineSegmentVector + lineSegment.point1;

            if (Utilities.AreIncident(linePlaneIntersectionPoint, lineSegment))
            {
                if (Utilities.AreIncident(linePlaneIntersectionPoint, triangle))
                    return Utilities.ShortestLineSegmentBetween(linePlaneIntersectionPoint, linePlaneIntersectionPoint);
                
                else
                    return Utilities.GetShortest(new Func<LineSegment3D>[]
                    {
                        () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge1),
                        () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge2),
                        () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge3)
                    });
            }

            else
                return Utilities.GetShortest(new Func<LineSegment3D>[]
                {
                    () => Utilities.ShortestLineSegmentBetween(lineSegment.point1, triangle),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment.point2, triangle),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge1),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge2),
                    () => Utilities.ShortestLineSegmentBetween(lineSegment, triangleEdge3)
                });
        }

        public static LineSegment3D ShortestLineSegmentBetween(Triangle3D triangle1, Triangle3D triangle2)
        {
            return Utilities.GetShortest(new Func<LineSegment3D>[]
            {
                () => Utilities.ShortestLineSegmentBetween(triangle1.Edge1, triangle2),
                () => Utilities.ShortestLineSegmentBetween(triangle1.Edge2, triangle2),
                () => Utilities.ShortestLineSegmentBetween(triangle1.Edge3, triangle2),
                () => Utilities.ShortestLineSegmentBetween(triangle2.Edge1, triangle1),
                () => Utilities.ShortestLineSegmentBetween(triangle2.Edge2, triangle1),
                () => Utilities.ShortestLineSegmentBetween(triangle2.Edge3, triangle1)
            });
        }
    }
}