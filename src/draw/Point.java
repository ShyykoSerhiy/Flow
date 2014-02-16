package draw;


/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:14
 */
public class Point {
    private double x;
    private double y;

    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public Point add(Point point) {
        return new Point(this.getX() + point.getX(), this.getY() + point.getY());
    }

    public Point minus(Point point) {
        return new Point(this.getX() - point.getX(), this.getY() - point.getY());
    }

    public Point divide(double scalar) {
        return new Point(this.getX() / scalar, this.getY() / scalar);
    }

    public Point multiply(double scalar) {
        return new Point(this.getX() * scalar, this.getY() * scalar);
    }

    public double multiply(Point point) {
        return this.getX() * point.getX() + this.getY() * point.getY();
    }
}
