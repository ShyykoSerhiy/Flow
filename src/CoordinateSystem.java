import processing.core.PApplet;

/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:05
 */
public class CoordinateSystem {
    private int width;
    private int height;
    private int maxX;
    private int maxY;
    private PApplet drawer;

    public CoordinateSystem(PApplet drawer,int width, int height, int maxX, int maxY) {
        this.width = width;
        this.height = height;
        this.maxX = maxX;
        this.maxY = maxY;
        this.drawer = drawer;
    }

    public void draw() {
        drawer.strokeWeight(1);
        drawer.fill(0, 0, 0);
        drawer.rect(0, 0, width, 5);
        drawer.text("X", width-20, 20);
        drawer.rect(0, 0, 5, height);
        drawer.text("Y", 10, height-10);

        int n = 30;
        float partX = maxX / ((float)n);
        float partY = maxY / ((float)n);
        int partW = (int) (width/((float)n));
        int partH = (int) (height/((float)n));
        for (int i = 0; i < n; i++){
            drawer.fill(0, 255, 0);
            drawer.rect(i*partW-5, 0, 5, 10);
            drawer.text(""+partX*i, i*partW-5, 30);
            drawer.rect(0, i*partH-5, 10, 5);
            drawer.text(""+partY*i, 30, i*partH-5);
        }
    }
}