using System.Windows;
using MathNet.Numerics;
using TinySpline;
namespace WpfApp1
{
    public class CatmullRom{
        public double DotProduct(Vector first, Vector other) => first.X * other.X + first.Y * other.Y;
        public double ToRadians(double angle) => (Math.PI / 180) * angle;
       
        public bool IsPointInCone(Point p1, Point p2, double angleDegrees, Point testPoint)
        {
            Vector direction = p2 - p1;
            Vector toTestPoint = testPoint - p2;

            double angleRadians = ToRadians(angleDegrees);
            double dotProduct = DotProduct(direction, toTestPoint);
            double angleToTestPoint = Math.Acos(dotProduct / (direction.Length * toTestPoint.Length));

            return angleToTestPoint <= angleRadians;
        }

        public Vector RotateVector(Vector v, double angle)
        {
            double cosTheta = Math.Cos(angle);
            double sinTheta = Math.Sin(angle);
            return new Vector(
                (v.X * cosTheta - v.Y * sinTheta),
                (v.X * sinTheta + v.Y * cosTheta)
            );
        }

        public List<Point> CRChain(List<Point> points, float alpha)
        {
            var chainedCoords = new List<Point>();
            for (int i = 0; i < points.Count - 3; i++)
            {
                chainedCoords.AddRange(CRPoint(points[i], points[i + 1], points[i + 2], points[i + 3], alpha));
            }
            return chainedCoords;
        }

        private double Distance(double ti, Point pi, Point pj, float alpha)
        {
            double xi = pi.X; double yi = pi.Y;
            double xj = pj.X; double yj = pj.Y;
            double dx = xj - xi; double dy = yj - yi;
            double l = Math.Sqrt(dx * dx + dy * dy);
            return ti + Math.Pow(l, alpha);
        }

        // Adapted wikipedia method https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline#Code_example_in_Python
        private List<Point> CRPoint(Point p0, Point p1, Point p2, Point p3, float alpha)
        {
            double t0 = 0;
            double t1 = Distance(t0, p0, p1, alpha);
            double t2 = Distance(t1, p1, p2, alpha);
            double t3 = Distance(t2, p2, p3, alpha);
            double[] t = Generate.LinearSpaced(1000, t1, t2);
            var points = new List<Point>();

            for (int i = 0; i < t.Length; i++)
            {
                var a1x = (t1 - t[i]) / (t1 - t0) * p0.X + (t[i] - t0) / (t1 - t0) * p1.X;
                var a1y = (t1 - t[i]) / (t1 - t0) * p0.Y + (t[i] - t0) / (t1 - t0) * p1.Y;

                var a2x = (t2 - t[i]) / (t2 - t1) * p1.X + (t[i] - t1) / (t2 - t1) * p2.X;
                var a2y = (t2 - t[i]) / (t2 - t1) * p1.Y + (t[i] - t1) / (t2 - t1) * p2.Y;

                var a3x = (t3 - t[i]) / (t3 - t2) * p2.X + (t[i] - t2) / (t3 - t2) * p3.X;
                var a3y = (t3 - t[i]) / (t3 - t2) * p2.Y + (t[i] - t2) / (t3 - t2) * p3.Y;

                var b1x = (t2 - t[i]) / (t2 - t0) * a1x + (t[i] - t0) / (t2 - t0) * a2x;
                var b1y = (t2 - t[i]) / (t2 - t0) * a1y + (t[i] - t0) / (t2 - t0) * a2y;

                var b2x = (t3 - t[i]) / (t3 - t1) * a2x + (t[i] - t1) / (t3 - t1) * a3x;
                var b2y = (t3 - t[i]) / (t3 - t1) * a2y + (t[i] - t1) / (t3 - t1) * a3y;

                var pointX = (t2 - t[i]) / (t2 - t1) * b1x + (t[i] - t1) / (t2 - t1) * b2x;
                var pointY = (t2 - t[i]) / (t2 - t1) * b1y + (t[i] - t1) / (t2 - t1) * b2y;
                points.Add(new Point(pointX, pointY));
            }
            return points;
        }

        public List<Point> CRLerp(List<Point> points, float alpha)
        {
            List<double> ps = new List<double>();

            foreach (var point in points)
            {
                ps.Add(point.X);
                ps.Add(point.Y);
            }
            var spline = BSpline.InterpolateCatmullRom(ps, 2, alpha);

            // Evaluate the spline at multiple points
            var resultPoints = new List<Point>();
            for (float u = 0.0f; u <= 1.0f; u += 1/(float)(100*(points.Count-1)))
            {
                var evaluated = spline.Eval(u).Result;
                resultPoints.Add(new Point(evaluated[0], evaluated[1]));
            }

            return resultPoints;
        }
        
        public bool IsSplineValid(List<Point> points, double angle)
        {
            for (int i = 2; i < points.Count; i++)
            {
                if (!IsPointInCone(points[i-2], points[i -1], angle, points[i]))
                {
                    return false;
                }
            }
            return true;
        }

    }
}