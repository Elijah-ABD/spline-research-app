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
        private List<Path> paths = new List<Path>();
        private CatmullRom cr = new CatmullRom();
        private float alpha = 0.5f;
        private bool draw = true;

        private UIElement selectedElement;
        private Point mouseOffset;

        public MainWindow() { InitializeComponent();}

        // Event Handlers
        public Dictionary<Path, Point> GetDict() => paths.Zip(coords, (k, v) => new { k, v }).ToDictionary(x => x.k, x => x.v);

        private void MouseClick(object sender, MouseEventArgs e)
        {
            if (draw)
            {
                var point = e.GetPosition(this);
                var path = new Path() { Fill = Brushes.LightBlue, Data = new EllipseGeometry(new Point(point.X, point.Y), 10, 10) };
                coords.Add(point);
                paths.Add(path);
                canvas.Children.Add(path);
                //Draw(point, el);
                if (coords.Count > 1) { SplineDrawing(Brushes.Orange); }
            }
        }
        private void sliderChange(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            alpha = (float)e.NewValue;
            if (coords.Count > 1) {SplineDrawing(Brushes.Orange);}
        }


        private void ClearButtonClick(object sender, RoutedEventArgs e) { clearCanvas(); coords.Clear(); paths.Clear(); }
        private void DrawButtonClick(object sender, RoutedEventArgs e) 
        { 
            draw = !draw;
            foreach (var path in paths)
            {
                EnableDrag(path);
            }
        }
        private void MapButtonClick(object sender, RoutedEventArgs e) => image.Visibility = image.Visibility
                                                                      == Visibility.Visible 
                                                                      ? Visibility.Hidden : Visibility.Visible;
        private void UndoButtonClick(object sender, RoutedEventArgs e)
        {
            switch (coords.Count)
            {
                case > 2:
                    coords.RemoveAt(coords.Count - 1);
                    SplineDrawing(Brushes.Orange);
                    break;

                case 2:
                    coords.RemoveAt(coords.Count - 1);
                    clearCanvas();
                    Draw(coords[0], 15, Brushes.LightGray);
                    break;

                default:
                    clearCanvas();
                    coords.Clear();
                    paths.Clear();
                    break;
            }
        }
        // Helper Functions
        private void Mirror(int index, Point p1, Point p2) => coords.Insert(index, new Point(2 * p1.X - p2.X, 2 * p1.Y - p2.Y));
        private void clearCanvas() => canvas.Children.Clear();

        // Drawing Functions
        private void SplineDrawing(Brush colour)
        {
            Mirror(0, coords[0], coords[1]);
            Mirror(coords.Count, coords[^1], coords[^2]);
            var splinePoints = cr.CRChain(coords, alpha);
            

            // Remove mirrored control points
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);

            clearCanvas();

            Draw(coords, 2, Brushes.LightGray);
            Draw(splinePoints, 5, colour);
            //foreach (var point in coords) { Draw(point, Brushes.LightGray); }
            foreach (var path in paths) { canvas.Children.Add(path); }

        }

        private void Draw(List<Point> points, int size, Brush colour)
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
    
        private void Draw(Point point, int size, Brush colour)
        {
            var el = new Ellipse
            {
                Width = size,
                Height = size,
                Fill = colour
            };
            Canvas.SetLeft(el, point.X-size/2);
            Canvas.SetTop(el, point.Y-size/2);
            canvas.Children.Add(el);
        }
        private void Draw(Point point, Brush colour)
        {
            Path path = new Path
            {
                Fill = colour,
                StrokeThickness = 1,
                Data = new EllipseGeometry(new Point(point.X, point.Y), 6, 6)
            };
            canvas.Children.Add(path);
        }

        private void Draw(Point start, Point end)
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

        public void DrawCone(Point p1, Point p2, double angleDegrees, double coneLength)
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



        private Nullable<Point> dragStart = null;

        private void Down(object sender, MouseEventArgs e) 
        { 
            var element = (UIElement)sender;
            dragStart = e.GetPosition(element);
            element.CaptureMouse();
        }
        private void Up(object sender, MouseEventArgs e) 
        { 
            var element = (UIElement)sender;
            dragStart = null;
            element.ReleaseMouseCapture();
        }

        private void Move(object sender, MouseEventArgs e) 
        { 
            if (dragStart != null && e.LeftButton == MouseButtonState.Pressed)
            {
                var element = (UIElement)sender;
                var p2 = e.GetPosition(this);
                Canvas.SetLeft(element, p2.X - dragStart.Value.X);
                Canvas.SetTop(element, p2.Y - dragStart.Value.Y);

                coords[paths.IndexOf((Path)element)] = p2;
                if (coords.Count > 1) SplineDrawing(Brushes.Orange);
                    
            }
        }

        private void EnableDrag(UIElement element) { 
                element.MouseDown += Down;
                element.MouseMove += Move;
                element.MouseUp += Up;
            }

        private void Mouse() 
        {

            var shapes = new UIElement[]
            {
                new Ellipse() { Fill = Brushes.DarkKhaki, Width = 100, Height = 100 },
                new Rectangle() { Fill = Brushes.LawnGreen, Width = 200, Height = 100 },
                new Path() { Fill = Brushes.LightBlue, Data = new EllipseGeometry(new Point(200,344),6,6)},
            };


        foreach (var shape in shapes) 
            {
                EnableDrag(shape);
                canvas.Children.Add(shape);
            }
        }
    }
}
