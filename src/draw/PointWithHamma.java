package draw;

/**
 * User: shyyko
 * Date: 16.02.14
 * Time: 12:42
 */
public class PointWithHamma extends Point {
    private double hamma;

    public PointWithHamma(Point point, double hamma) {
        super(point.getX(), point.getY());
        this.hamma = hamma;
    }

    public PointWithHamma(double x, double y, double hamma) {
        super(x, y);
        this.hamma = hamma;
    }

    public double getHamma() {
        return hamma;
    }

    public void setHamma(double hamma) {
        this.hamma = hamma;
    }


}
