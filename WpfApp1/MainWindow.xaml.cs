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

        private List<(double x, double y)> coords = new List<(double x, double y)>();
        private List<(double x, double y)> splineCoords = new List<(double x, double y)>();
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
            coords.Add((point.X, point.Y));

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
            for (double i = 0; i < coords.Count; i+=0.001)
            {
                //splineCoords.Add(coords[i]);
 
                    GetSplineCoord(i);
            }

        }
        private void GetSplineCoord(double t, bool loop=false)
        {
            int p0, p1, p2, p3;

            if (!loop)
            {
                p1 = (int)t + 1;
                p2 = p1 + 1;
                p3 = p2 + 1;
                p0 = p1 - 1;
            }
            else
            {
                p1 = (int)t;
                p2 = (p1 + 1) % coords.Count;
                p3 = (p2 + 1) % coords.Count;
                p0 = p1 > 1 ? p1 - 1 : coords.Count - 1;
            }

            t = t - (int)t;

            double tSq = t * t;
            double tCu = tSq * t;

            double q0, q1, q2, q3;

            q0 = -tCu + 2 * tSq - t;
            q1 = 3 * tCu - 5 * tSq + 2;
            q2 = -3 * tCu + 4 * tSq + t;
            q3 = tCu - tSq;

            double tX = 0.5 * (coords[p0].Item1 * q0 + coords[p1].Item1 * q1 + coords[p2].Item1 * q2 + coords[p3].Item1 * q3);
            double tY = 0.5 * (coords[p0].Item2 * q0 + coords[p1].Item2 * q1 + coords[p2].Item2 * q2 + coords[p3].Item2 * q3);

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

        private (double,double) getDouble(string s){
            var strings = s.Split(',');
            double x = Convert.ToDouble(strings[0]);
            double y = Convert.ToDouble(strings[1]);
            return (x, y);
        } 

        private void clearCanvas(){
            while (canvas.Children.Count > 2){
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
                    coords.Add(getDouble(tb.Text));
                }
            }
            clearCanvas();
            DrawSpline();
            coords.Clear();
            }
        }
    }
}
