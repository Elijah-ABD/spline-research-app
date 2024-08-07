using System.Security.Cryptography.X509Certificates;
using System.Windows;
using System.Windows.Automation;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using MathNet.Numerics;
using Microsoft.VisualBasic;

namespace WpfApp1
{
    struct Coord
    {
        public Coord(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
        public double x, y;
    }
    public partial class MainWindow : System.Windows.Window
    {
        private List<Coord> coords = new List<Coord>();

        public MainWindow()
        {
            InitializeComponent();
            canvas.Children.Add(new Line
            {
                X1 = 288,
                Y1 = 288,
                X2 = 288,
                Y2 = 86,
                Stroke = Brushes.White,
                StrokeThickness = 2
            });
            coords.Add(new Coord(288, 86));
            coords.Add(new Coord(288, 33));
            //coords.Add(new Coord(786, 37));
        }

        // Event Handlers

        private void MouseClickHandler(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            var newPoint = new Coord(point.X, point.Y);
            coords.Add(newPoint);

            //DrawPoints(new List<Coord>{newPoint}, 10, Brushes.Black);
            this.canvas.Children.Add(
                   new Path()
                   {
                       Stroke = Brushes.LightGray,
                       StrokeThickness = 0,
                       Data = new EllipseGeometry(point, 6, 6)
                   }
               );

            if (coords.Count > 1)
            {
                Mirror();
                clearCanvas();
                //SplineDrawing(Brushes.Red, 0);
                //SplineDrawing(Brushes.Blue, 1);

                //DrawPoints(Uniform(), 2, Brushes.Red, true);
                SplineDrawing(Brushes.Yellow, 0.5f);
                double angleOfView = 60;
                double distance = 500;
                // Calculate the angle of view in radians
                double angleOfViewRadians = angleOfView * Math.PI / 180;
                // Calculate the points of the cone
                var centerPoint = new Point(coords[^1].x, coords[^1].y);
                var point1 = new Point(centerPoint.X + distance * Math.Cos(angleOfViewRadians / 2),
                centerPoint.Y + distance * Math.Sin(angleOfViewRadians / 2));
                var point2 = new Point(centerPoint.X + distance * Math.Cos(-angleOfViewRadians / 2),
                centerPoint.Y + distance * Math.Sin(-angleOfViewRadians / 2));
                // Add the cone to the canvas
                canvas.Children.Add(new Polygon
                {
                Points = new PointCollection { centerPoint, point1, point2 },
                Stroke = Brushes.White,
                StrokeThickness = 0
                });


            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            var button = (Button)sender;
            switch (button.Name)
            {
                case "LoopButton":
                    if (coords.Count > 3)
                    {
                        DrawPoints(new List<Coord> { new(coords[^1].x, coords[^1].y), new(coords[0].x, coords[0].y) }, 1, Brushes.LightGray);
                        DrawPoints(CRPoint(coords[1], coords[0], coords[^1], coords[^2], 0.5f), 2, Brushes.Yellow);

                    }
                    break;

                case "ClearButton":
                    clearCanvas();
                    while (coords.Count > 2)
                    {
                        coords.RemoveAt(coords.Count - 1);
                    }
                    break;
            }


        }

        // Helper Functions
        private void Mirror()
        {
            var p1 = coords[0];
            var p2 = coords[1];

            var dx = p2.x - p1.x;
            var dy = p2.y - p1.y;

            var first = new Coord(p1.x - dx, p1.y - dy);
            coords.Insert(0, first);

            p1 = coords[^2];
            p2 = coords[^1];

            dx = p2.x - p1.x;
            dy = p2.y - p1.y;

            var last = new Coord(p2.x + dx, p2.y + dy);
            coords.Add(last);

        }


        private (double, double) getDouble(string s)
        {
            var strings = s.Split(',');
            double x = Convert.ToDouble(strings[0]);
            double y = Convert.ToDouble(strings[1]);
            return (x, y);
        }

        private void clearCanvas()
        {
            while (canvas.Children.Count > 1)
            {
                canvas.Children.RemoveAt(canvas.Children.Count - 1);
            } 

        }

        // TODO: ADD FILL FUNCTIONALITY
        private void SplineDrawing(Brush colour, float alpha = 0)
        {
            var splinePoints = CRChain(coords, alpha);
            //var zero = coords[0];
            //var end = coords[^1];
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);
            DrawPoints(coords, 2, Brushes.LightGray);
            //DrawPoints(coords, 10, Brushes.Black);
            foreach (var point in coords)
            {
                this.canvas.Children.Add(
                    new Path()
                    {
                        Stroke = Brushes.LightGray,
                        StrokeThickness =0,
                        Data = new EllipseGeometry((new Point(point.x, point.y)), 6, 6)
                    }
                );
            }
            DrawPoints(splinePoints, 2, colour);
            //coords.Insert(0,zero);
            //coords.Add(end);

        }

        private List<Coord> CRChain(List<Coord> points, float alpha)
        {
            var chainedCoords = new List<Coord>();
            for (int i = 0; i < points.Count - 3; i++)
            {
                chainedCoords.AddRange(CRPoint(points[i], points[i + 1], points[i + 2], points[i + 3], alpha));
            }
            return chainedCoords;
        }

        private double Distance(double ti, Coord pi, Coord pj, float alpha)
        {
            double xi = pi.x; double yi = pi.y;
            double xj = pj.x; double yj = pj.y;
            double dx = xj - xi; double dy = yj - yi;
            double l = Math.Sqrt(dx * dx + dy * dy);
            return ti + Math.Pow(l, alpha);
        }

        // Adapted wikipedia method https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline#Code_example_in_Python
        private List<Coord> CRPoint(Coord p0, Coord p1, Coord p2, Coord p3, float alpha)
        {
            double t0 = 0;
            double t1 = Distance(t0, p0, p1, alpha);
            double t2 = Distance(t1, p1, p2, alpha);
            double t3 = Distance(t2, p2, p3, alpha);
            double[] t = Generate.LinearSpaced(1000, t1, t2);
            var points = new List<Coord>();
            for (int i = 0; i < t.Length; i++)
            {
                var a1x = (t1 - t[i]) / (t1 - t0) * p0.x + (t[i] - t0) / (t1 - t0) * p1.x;
                var a1y = (t1 - t[i]) / (t1 - t0) * p0.y + (t[i] - t0) / (t1 - t0) * p1.y;

                var a2x = (t2 - t[i]) / (t2 - t1) * p1.x + (t[i] - t1) / (t2 - t1) * p2.x;
                var a2y = (t2 - t[i]) / (t2 - t1) * p1.y + (t[i] - t1) / (t2 - t1) * p2.y;

                var a3x = (t3 - t[i]) / (t3 - t2) * p2.x + (t[i] - t2) / (t3 - t2) * p3.x;
                var a3y = (t3 - t[i]) / (t3 - t2) * p2.y + (t[i] - t2) / (t3 - t2) * p3.y;

                var b1x = (t2 - t[i]) / (t2 - t0) * a1x + (t[i] - t0) / (t2 - t0) * a2x;
                var b1y = (t2 - t[i]) / (t2 - t0) * a1y + (t[i] - t0) / (t2 - t0) * a2y;

                var b2x = (t3 - t[i]) / (t3 - t1) * a2x + (t[i] - t1) / (t3 - t1) * a3x;
                var b2y = (t3 - t[i]) / (t3 - t1) * a2y + (t[i] - t1) / (t3 - t1) * a3y;

                var pointX = (t2 - t[i]) / (t2 - t1) * b1x + (t[i] - t1) / (t2 - t1) * b2x;
                var pointY = (t2 - t[i]) / (t2 - t1) * b1y + (t[i] - t1) / (t2 - t1) * b2y;
                points.Add(new Coord { x = pointX, y = pointY });
            }
            return points;
        }

        private void DrawPoints(List<Coord> points, int size, Brush colour, bool tracing = false)
        {
            if (!tracing)
            {
                foreach (var point in points)
                {
                    Ellipse marker = new Ellipse
                    {
                        Height = size,
                        Width = size,
                        Fill = colour,
                    };
                    Canvas.SetLeft(marker, point.x - marker.Width / 2);
                    Canvas.SetTop(marker, point.y - marker.Height / 2);
                    canvas.Children.Add(marker);
                }
            }
            else
            {
                for (int i = 0; i < points.Count - 1; i++)
                {
                    Line line = new Line
                    {
                        X1 = points[i].x,
                        Y1 = points[i].y,
                        X2 = points[i + 1].x,
                        Y2 = points[i + 1].y,
                        Stroke = colour,
                        StrokeThickness = size
                    };
                    canvas.Children.Add(line);
                }
            }

        }
        // old
                private void DrawSpline()
        {
            var splineCoords = new List<Coord>();
            for (int i = 0; i < coords.Count - 3; i++)
            {
                splineCoords.Add(coords[i]);
                var points = new List<Coord>() { coords[i], coords[i + 1], coords[i + 2], coords[i + 3] };

                for (double t = 0; t < 1; t += 0.001)
                {
                    GetSplineCoord(t, points);
                }
            }


        }
        private void GetSplineCoord(double t, List<Coord> points)
        {
            int p0, p1, p2, p3;

            p0 = 0;
            p1 = 1;
            p2 = 2;
            p3 = 3;

            double tSq = t * t;
            double tCu = tSq * t;

            double q0, q1, q2, q3, tX, tY;

            q0 = -tCu + 3*tSq - 3*t +1;
            q1 = 3*tCu - 6*tSq + 4;
            q2 = -3*tCu + 3*tSq + 3*t + 1;
            q3 = tCu;

            tX = (1/6.0) * (points[p0].x * q0 + points[p1].x * q1 + points[p2].x * q2 + points[p3].x * q3);
            tY = (1/6.0) * (points[p0].y * q0 + points[p1].y * q1 + points[p2].y * q2 + points[p3].y * q3);

            Ellipse marker = new Ellipse
            {
                Height = 2,
                Width = 2,
                Fill = Brushes.Blue,
            };
            Canvas.SetLeft(marker, tX - marker.Width / 2);
            Canvas.SetTop(marker, tY - marker.Height / 2);
            canvas.Children.Add(marker);
        }
        // malu
        private List<Coord> Uniform()
        {
            var controlPoints2 = coords;
            controlPoints2.Insert(0, new Coord(2 * coords[0].x - coords[1].x, 2 * coords[0].y - coords[1].y));
            controlPoints2.Insert(controlPoints2.Count, new Coord(2 * coords[coords.Count - 1].x - coords[coords.Count - 2].x, 2 * coords[coords.Count - 1].y - coords[coords.Count - 2].y));

            double ts = 3.0; //tension
            MatrixNxM cubicBlendMatrix = (1.0 / 6.0) * new MatrixNxM(
                new double[4, 4]
                {
                    {2-ts,  6-ts, ts-6,  ts-2},
                    {2*ts-3,ts-9, 9-2*ts,3-ts},
                    {-ts,   0,    ts,    0},
                    { 1,    4,    1,     0}
                });

            int Segments = controlPoints2.Count - 4 + 1;

            List<Coord> curve = new List<Coord>();

            for (int s = 0; s < Segments; ++s)
            {
                for (double t = 0; t < 1.001; t += 0.001)
                {
                    MatrixNxM tvector = new MatrixNxM(new double[1, 4] { { t * t * t, t * t, t, 1 } });
                    MatrixNxM basis = tvector * cubicBlendMatrix;

                    Coord p = new Coord();
                    p.x = 0; p.y = 0;
                    for (int i = 0; i < 4; ++i)
                    {
                        p.x += basis[0, i] * controlPoints2[s + i].x;
                        p.y += basis[0, i] * controlPoints2[s + i].y;
                    }
                    curve.Add(p);
                }
            }

            return curve;
        }
    }
}
