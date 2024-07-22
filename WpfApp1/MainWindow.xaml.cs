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
using System.Windows.Navigation;
using System.Windows.Shapes;


namespace WpfApp1
{
    public partial class MainWindow : Window
    {
        private List<(double x, double y)> coords = new List<(double x, double y)> ();
        private List<(double x, double y)> splineCoords = new List<(double x, double y)> ();
        public MainWindow()
        {
            InitializeComponent();
        }


        private void MouseMoveHandler(object sender, MouseEventArgs e)
        {
            textBox1.Text = $"{e.GetPosition(this)}";
        }

        private void MouseClickHandler(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            coords.Add((point.X,point.Y));
            textBox2.Text = $"{String.Join("; ", coords)}";

            Ellipse marker = new Ellipse
            {
                Height = 10,
                Width = 10,
                Fill = Brushes.Red,
            };
            Canvas.SetLeft(marker, point.X - marker.Width / 2);
            Canvas.SetTop(marker, point.Y - marker.Height / 2);
            canvas.Children.Add(marker);

            
        }

        private void OnKeyDownHandler(object sender, KeyEventArgs e)
        {
            if ((e.Key == Key.Return) || (e.Key == Key.Space))
            {
                DrawSplin();
            }

            if (e.Key == Key.R) 
            {
                canvas.Children.Clear();
                coords.Clear();
            }

            if (e.Key == Key.Escape)
            {
                Close();
            }
        }

        private void DrawSpline()
        {
            List<Line> lines = new List<Line>();
            if (coords.Count < 4)
            {
                return;
            }
            else
            {
                for (int i = 0; i < coords.Count - 1; i++)
                {
                    Line line = new Line
                    {
                        X1 = coords[i].Item1,
                        Y1 = coords[i].Item2,
                        X2 = coords[i+1].Item1,
                        Y2 = coords[i+1].Item2,
                        Stroke = Brushes.Black,
                        StrokeThickness = 4,
                    };
                    lines.Add(line);
                }
                foreach (Line line in lines)
                {

                    canvas.Children.Add(line);
                }
            }

         

        }

        private void DrawSplin()
        {
            for (int i = 0; i < coords.Count - 3; i++)
            {
                splineCoords.Add(coords[i]);
                var points = new List<(double, double)>() { coords[i], coords[i+1], coords[i+2], coords[i+3] };

                for (double t = 0; t < 1; t += 0.005)
                {
                    GetSplineCoord(t, points);
                }
            }

        }
        private void  GetSplineCoord(double t, List<(double, double)> points )
        {
            int p0, p1, p2, p3;
            
            p0 = 0;
            p1 = 1;
            p2 = 2;
            p3 = 3;

            double tSq = t * t;
            double tCu = tSq * t;

            double q0, q1, q2, q3;

            q0 = -tCu + 2*tSq - t; 
            q1 = 3*tCu - 5*tSq +2; 
            q2 = -3*tCu + 4*tSq + t; 
            q3 = tCu - tSq; 

            double tX = 0.5 * (points[p0].Item1 * q0 + points[p1].Item1 * q1 + points[p2].Item1 * q2 + points[p3].Item1 * q3);
            double tY = 0.5 * (points[p0].Item2 * q0 + points[p1].Item2 * q1 + points[p2].Item2 * q2 + points[p3].Item2 * q3);

            Ellipse marker = new Ellipse
            {
                Height = 10,
                Width = 10,
                Fill = Brushes.Black,
            };
            Canvas.SetLeft(marker, tX - marker.Width / 2);
            Canvas.SetTop(marker, tY - marker.Height / 2);
            canvas.Children.Add(marker);
        }

    }
}