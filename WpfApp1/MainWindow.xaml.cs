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


namespace WpfApp1
{
    public partial class MainWindow : System.Windows.Window
    {

        private List<Coord> coords = new List<Coord>();
        private List<Coord> splineCoords = new List<Coord>();
        struct Coord
        {
            public double x, y;
        }
        public MainWindow()
        {
            InitializeComponent();
        }


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


        }


        private void OnKeyDownHandler(object sender, KeyEventArgs e)
        {
            if ((e.Key == Key.Return) || (e.Key == Key.Space))
            {
                DrawSpline();
                SplineDrawing();
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
        //private void DrawSpline()
        //{
        //    for (int i = 0; i < coords.Count - 3; i++)
        //    {
        //        splineCoords.Add(coords[i]);
        //        var points = new List<(double, double)>() { coords[i], coords[i + 1], coords[i + 2], coords[i + 3] };

        //        for (double t = 0; t < 1; t += 0.005)
        //        {
        //            GetSplineCoord(t, points, true);
        //        }
        //    }

        //}


        private void DrawSpline()
        {
            for (int i = 0; i < coords.Count - 3; i++)
            {
                //splineCoords.Add(coords[i]);
                for (double t = 0; t < 1; t += 0.001)
                {
                    var point = GetSplineCoord(t, coords[i], coords[i + 1], coords[i + 2], coords[i + 3]);
                    Ellipse marker = new Ellipse
                    {
                        Height = 2,
                        Width = 2,
                        Fill = Brushes.Blue,
                    };
                    Canvas.SetLeft(marker, point.x - marker.Width / 2);
                    Canvas.SetTop(marker, point.y - marker.Height / 2);
                    canvas.Children.Add(marker);
                }
            }

        }
        private Coord GetSplineCoord(double t, Coord p0, Coord p1, Coord p2, Coord p3, bool loop = false)
        {
            var point = new Coord();
            var t2 = t * t;
            var t3 = t2 * t;


            point.x = 0.5f * ((2.0f * p1.x)
                + (-p0.x + p2.x) * t
                + (2.0f * p0.x - 5.0f * p1.x
                + 4 * p2.x - p3.x) * t2 +
                (-p0.x + 3.0f * p1.x - 3.0f * p2.x + p3.x) * t3);

            point.y = 0.5f * ((2.0f * p1.y)
                + (-p0.y + p2.y) * t
                + (2.0f * p0.y - 5.0f * p1.y + 4 * p2.y - p3.y) * t2
                + (-p0.y + 3.0f * p1.y - 3.0f * p2.y + p3.y) * t3);


            return point;

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
            while (canvas.Children.Count > 2)
            {
                canvas.Children.RemoveAt(canvas.Children.Count - 1);

            }
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
                DrawSpline();
                coords.Clear();
            }
        }
        private void SplineDrawing()
        {
            var centripetalPoints = CRChain(coords, 1);
            DrawPoints(coords, 10, Brushes.Black);
            DrawPoints(centripetalPoints, 2, Brushes.Red);
            
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

    }
}
