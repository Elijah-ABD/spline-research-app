using Accessibility;
using System.Drawing;
using System.Security.Cryptography.Pkcs;
using System.Text;
using System.Transactions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MathNet.Numerics;
using System.Security.Cryptography.Xml;
using System.Security.Policy;


namespace WpfApp1
{
    struct Coord
    {
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

        private void MouseMoveHandler(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            textBox1.Text = $"x: {Math.Round(point.X, 2)}, y: {Math.Round(point.Y, 2)}";
        }

        private void MouseClickHandler(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            Coord p;
            coords.Add(new Coord { x = point.X, y = point.Y });

            Ellipse marker = new Ellipse
            {
                Height = 10,
                Width = 10,
                Fill = Brushes.Black,
            };
            Canvas.SetLeft(marker, point.X - marker.Width / 2);
            Canvas.SetTop(marker, point.Y - marker.Height / 2);
            canvas.Children.Add(marker);

            if (coords.Count > 1)
            {
                Mirror();
                clearCanvas();
                SplineDrawing(Brushes.Red, 0.5f);
            }


        }


        private void OnKeyDownHandler(object sender, KeyEventArgs e)
        {
            if ((e.Key == Key.Return) || (e.Key == Key.Space))
            {
                if (coords.Count > 3)
                {
                    DrawPoints(CRPoint(coords[1], coords[0], coords[^1], coords[^2], 0.5f), 2, Brushes.Red);
                }
               
            }

            if (e.Key == Key.R)
            {
                clearCanvas();
                coords.Clear();
            }

            if (e.Key == Key.Escape)
            {
                Close();
            }
        }

        // Helper Functions
        private void Mirror()
        {
                var p1 = coords[0];
                var p2 = coords[1];

                var dx = p2.x - p1.x;
                var dy = p2.y - p1.y;

                var first = new Coord { x = p1.x-dx, y = p1.y-dy };
                coords.Insert(0, first);

                p1 = coords[^2];
                p2 = coords[^1];
                
                dx = p2.x - p1.x;
                dy = p2.y - p1.y;
                
                var last = new Coord { x = p2.x + dx, y = p2.y + dy };
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

        private void TextBox_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Enter)
            {
                foreach (UIElement element in sp1.Children)
                {
                    if (element is TextBox tb)
                    {
                        coords.Add(new Coord { x = getDouble(tb.Text).Item1, y = getDouble(tb.Text).Item2 });
                    }
                }
                clearCanvas();
                SplineDrawing(Brushes.Blue);
                coords.Clear();
            }
        }

        // Spline draw functionality
        private void SplineDrawing(Brush colour, float alpha=0, bool fill=false)
        {
            var centripetalPoints = CRChain(coords, alpha, fill);
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);
            DrawPoints(coords, 10, Brushes.Black);
            DrawPoints(centripetalPoints, 2, colour);
            
        }

        private List<Coord> CRChain(List<Coord> points, float alpha, bool fill=false)
        {
            var chainedCoords = new List<Coord>();
            if (fill)
            {
                for (int i = 0; i < points.Count; i++)
                {
                    chainedCoords.AddRange(CRPoint(points[i], points[(i + 1) % points.Count], points[(i + 2) % points.Count], points[(i + 3) % points.Count], alpha));
                }
            }
            else
            {
                for (int i = 0; i < points.Count - 3; i++)
                {
                    chainedCoords.AddRange(CRPoint(points[i], points[i + 1], points[i + 2], points[i + 3], alpha));
                }
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

        private void DrawPoints(List<Coord> points,int size, Brush colour) 
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

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            if (coords.Count > 3)
            {
                DrawPoints(CRPoint(coords[1], coords[0], coords[^1], coords[^2], 0.5f), 2, Brushes.Red);
            }
        }
    }
}
