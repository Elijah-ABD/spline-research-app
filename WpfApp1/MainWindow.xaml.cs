using System.Diagnostics;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
namespace WpfApp1
{
    public partial class MainWindow : System.Windows.Window
    {
        private List<Point> coords = new List<Point>();
        private CatmullRom cr = new CatmullRom();
        private float alpha = 0;
        private double magnitude = 0.5;
        private Point first;

        public MainWindow()
        {
            InitializeComponent();

        }

        // Event Handlers
        private void MouseClick(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            coords.Add(point);
            Draw(point, 6, Brushes.LightGray);
            Connect();
        }

        private void Connect()
        {
            if (coords.Count == 2)
            {
                var line = new Line
                {
                    X1 = coords[0].X,
                    X2 = coords[1].X,
                    Y1 = coords[0].Y,
                    Y2 = coords[1].Y,
                    StrokeThickness = 2,
                    Stroke = Brushes.Bisque,
                };

                canvas.Children.Add(line);
                var p = new Point(line.X1 + magnitude * (line.X2 - line.X1), line.Y1 + magnitude * (line.Y2 - line.Y1));

                Draw(p, 6, Brushes.LightGreen);
                first = coords[0];
                coords.RemoveAt(0);

                coords.Insert(0, p);
            }
            if (coords.Count >= 3)
            {

                coords.RemoveAt(0);
                coords.Insert(0, first);


                var line = new Line
                {
                    X1 = coords[0].X,
                    X2 = coords[1].X,
                    Y1 = coords[0].Y,
                    Y2 = coords[1].Y,
                    StrokeThickness = 2,
                    Stroke = Brushes.Bisque,
                };
                clearCanvas();
                canvas.Children.Add(line);
                var p = new Point(line.X1 + magnitude * (line.X2 - line.X1), line.Y1 + magnitude * (line.Y2 - line.Y1));
                Draw(first, 6, Brushes.LightGray);
                Draw(p, 6, Brushes.LightGreen);
                coords.RemoveAt(0);
                coords.Insert(0, p);
                SplineDrawing(Brushes.Orange);
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            switch (((Button)sender).Content)
            {
                case "Clear":
                    clearCanvas();
                    coords.Clear();
                    break;

                case "Toggle Map":
                    image.Visibility = image.Visibility == Visibility.Visible ? Visibility.Hidden : Visibility.Visible;
                    break;
            }
        }

        // Helper Functions
        private void MirrorControlPoints()
        {
            // Mirroring the first and last control point
            var p1 = coords[0];
            var p2 = coords[1];
            //coords.Insert(0, new Point(2*p1.X-p2.X, 2* p1.Y - p2.Y));

            p1 = coords[^2];
            p2 = coords[^1];

            var last = new Point(2 * p2.X - p1.X, 2 * p2.Y - p1.Y);
            coords.Add(last);
        }

        private void clearCanvas() => canvas.Children.Clear();

        private void SplineDrawing(Brush colour)
        {
            MirrorControlPoints();
            var splinePoints = cr.CRChain(coords, alpha);

            // Remove mirrored control points
            //coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);

            //clearCanvas();
            //Draw(coords, 2, Brushes.LightGray);
            foreach (var point in coords)
            {
                Draw(point, 6, Brushes.LightGray);
            }
            Draw(splinePoints, 2, colour);

        }

        private void Draw(List<Point> points, int size, Brush colour)

        {
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

        }

        private void Draw(Point point, int size, Brush colour)
        {
            canvas.Children.Add(
                new Path()
                {
                    Stroke = colour,
                    StrokeThickness = 1,
                    Data = new EllipseGeometry(new Point(point.X, point.Y), size, size)
                }
            );
        }

        private void sliderChange(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            //alpha = (float)e.NewValue;
            //if (coords.Count > 1){SplineDrawing(Brushes.Orange);}
            magnitude = e.NewValue;
            coords.RemoveAt(0);
            coords.Insert(0, first);


            var line = new Line
            {
                X1 = coords[0].X,
                X2 = coords[1].X,
                Y1 = coords[0].Y,
                Y2 = coords[1].Y,
                StrokeThickness = 2,
                Stroke = Brushes.Bisque,
            };
            clearCanvas();
            canvas.Children.Add(line);
            var p = new Point(line.X1 + magnitude * (line.X2 - line.X1), line.Y1 + magnitude * (line.Y2 - line.Y1));
            Draw(first, 6, Brushes.LightGray);
            Draw(p, 6, Brushes.LightGreen);
            coords.RemoveAt(0);
            coords.Insert(0, p);
            SplineDrawing(Brushes.Orange);

        }

    }
}
