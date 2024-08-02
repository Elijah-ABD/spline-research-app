using System.Security.Cryptography.X509Certificates;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using MathNet.Numerics;

namespace WpfApp1
{
    struct Coord
    {
        public Coord(double x, double y){
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
        }

        // Event Handlers

        private void MouseClickHandler(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            var newPoint = new Coord(point.X, point.Y);
            coords.Add(newPoint);

            DrawPoints(new List<Coord>{newPoint}, 10, Brushes.Black);

            if (coords.Count > 1)
            {
                Mirror();
                clearCanvas();
                SplineDrawing(0.5f);

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
                            DrawPoints(new List<Coord>{new(coords[^1].x, coords[^1].y), new(coords[0].x, coords[0].y)}, 2, Brushes.LightGray, true);
                            DrawPoints(CRPoint(coords[1], coords[0], coords[^1], coords[^2], 0.5f), 2, Brushes.Red);
                            
                        }
                    break;
                
                case "ClearButton":
                    this.canvas.Children.Clear();
                    coords.Clear(); 
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

                var first = new Coord(p1.x-dx, p1.y-dy);
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
            this.canvas.Children.Clear();
        }

        // TODO: ADD FILL FUNCTIONALITY
        private void SplineDrawing(float alpha=0)
        {
            var splinePoints = CRChain(coords, alpha);
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);
            DrawPoints(coords, 2, Brushes.LightGray, true);
            DrawPoints(coords, 10, Brushes.Black);
            DrawPoints(splinePoints, 2, Brushes.Red);
            
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

        private void DrawPoints(List<Coord> points,int size, Brush colour, bool tracing=false) 
        {
            if (!tracing){
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
                        X1 = points[i].x, Y1 = points[i].y,
                        X2 = points[i+1].x, Y2 = points[i+1].y,
                        Stroke = colour,
                        StrokeThickness=size
                    };
                    canvas.Children.Add(line);
                }
            }

        }
    }
}
