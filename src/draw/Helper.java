package draw;

/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:13
 */
public class Helper {
    private int width;
    private int height;
    private float maxX;
    private float maxY;

    public Helper(int width, int height, float maxX, float maxY) {
        this.width = width;
        this.height = height;
        this.maxX = maxX;
        this.maxY = maxY;
    }

    public int getWidth() {
        return width;
    }

    public void setWidth(int width) {
        this.width = width;
    }

    public int getHeight() {
        return height;
    }

    public void setHeight(int height) {
        this.height = height;
    }

    public float getMaxX() {
        return maxX;
    }

    public void setMaxX(float maxX) {
        this.maxX = maxX;
    }

    public float getMaxY() {
        return maxY;
    }

    public void setMaxY(float maxY) {
        this.maxY = maxY;
    }

    public Point getCoordinatePoint(Point drawPoint) {
        return new Point(drawPoint.getX() / width * maxX, drawPoint.getY() / height * maxY);
    }

    public Point getDrawPoint(Point coordinatePoint) {
        return new Point(coordinatePoint.getX() / maxX * width, coordinatePoint.getY() / maxY * height);
    }
}
