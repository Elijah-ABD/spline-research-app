using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Shapes;
using System.Windows;

namespace WpfApp1
{
    public class Drawer
    {
        private readonly CatmullRom cr = new CatmullRom();
        private readonly Canvas canvas;
        private readonly List<Point> coords;
        private readonly List<Path> paths;

        public Drawer(Canvas canvas, List<Point> coords, List<Path> paths)
        {
            this.canvas = canvas;
            this.coords = coords;
            this.paths = paths;
        }
        public void clearCanvas() => canvas.Children.Clear();
        public void SplineDrawing(Brush colour, float alpha, bool draw, int splineType=1)
        {
            if (splineType == 3 && coords.Count < 4) return; 
            var splinePoints = cr.CRLerp(coords, alpha, splineType);
            var polyline = new Polyline{Stroke = colour, StrokeThickness = 5};
            clearCanvas();

            Draw(coords, 2, Brushes.LightGray);
            foreach (var point in splinePoints) { polyline.Points.Add(point); }
            canvas.Children.Add(polyline);
            foreach (var path in paths) { canvas.Children.Add(path); }
            
            if (draw && coords.Count > 1) Draw(coords[^2], coords[^1], 45, 2000); // Arbitrary Values
        }

        public void Draw(List<Point> points, int size, Brush colour)
        {
            for (int i = 0; i < points.Count - 1; i++)
            {
                Line line = new Line
                {
                    X1 = points[i].X,
                    Y1 = points[i].Y,
                    X2 = points[i + 1].X,
                    Y2 = points[i + 1].Y,
                    Stroke = colour,
                    StrokeThickness = size
                };
                canvas.Children.Add(line);
            }
        }

        public Path Draw(Point point, int size, Brush colour)
        {
            Path path = new Path
            {
                Fill = colour,
                StrokeThickness = 1,
                Data = new EllipseGeometry(new Point(point.X, point.Y), size, size)
            };
            canvas.Children.Add(path);
            return path;
        }

        public void Draw(Point start, Point end)
        {
            Line line = new Line
            {
                X1 = start.X,
                Y1 = start.Y,
                X2 = end.X,
                Y2 = end.Y,
                Stroke = Brushes.LightGray,
                StrokeThickness = 2
            };
            canvas.Children.Add(line);
        }

        public void Draw(Point p1, Point p2, double angleDegrees, double coneLength)
        {
            Vector direction = p2 - p1;
            direction.Normalize();
            double angleRadians = cr.ToRadians(angleDegrees);

            Vector scaledDirection = direction * coneLength;
            Vector leftEdge = cr.RotateVector(scaledDirection, -angleRadians);
            Vector rightEdge = cr.RotateVector(scaledDirection, angleRadians);

            Point leftPoint = p2 + leftEdge;
            Point rightPoint = p2 + rightEdge;

            Draw(p2, leftPoint);
            Draw(p2, rightPoint);
        }
    }
}
