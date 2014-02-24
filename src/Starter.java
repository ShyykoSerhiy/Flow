import draw.Helper;
import processing.core.PApplet;
import processing.core.PVector;
import shape.ShShape;
import shape.Shape;
import shape.ShapeDrawer;
import solver.Solver;

/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:03
 */
public class Starter extends PApplet {
    private static final float COORDINATE_MULTIPLIER = 1f;
    float zoom;
    // A vector to store the offset from the center
    PVector offset;
    // The previous offset
    PVector poffset;
    // A vector for the mouse position
    PVector mouse;
    CoordinateSystem coordinateSystem;
    Shape shape;
    ShapeDrawer shapeDrawer;
    Solver solver;
    final int width = 640;
    final int height = 640;
    boolean toDrawPhi = true;
    boolean toDrawVectors = true;
    // Figure figure;
    // Solver solver;

    @Override
    public void setup() {
        size(width, height);
        zoom = 1.0f;
        offset = new PVector(0, 0);
        poffset = new PVector(0, 0);
        coordinateSystem = new CoordinateSystem(this, (int) (width * COORDINATE_MULTIPLIER), (int) (height * COORDINATE_MULTIPLIER),
                (int) (3 * COORDINATE_MULTIPLIER), (int) (3 * COORDINATE_MULTIPLIER));
        shape = new ShShape(204);
        solver = new Solver(shape);
        shapeDrawer = new ShapeDrawer(new Helper((int) (width * COORDINATE_MULTIPLIER), (int) (height * COORDINATE_MULTIPLIER),
                (int) (3 * COORDINATE_MULTIPLIER), (int) (3 * COORDINATE_MULTIPLIER)), shape, this, solver);
        smooth();
    }

    public void draw() {
        background(255);
        pushMatrix();
        // Everything must be drawn relative to center
        translate(0, 0);

        // Use scale for 2D "zoom"
        scale(zoom);
        // The offset (note how we scale according to the zoom)
        translate(offset.x / zoom, offset.y / zoom);

        // An arbitrary design so that we have something to see!
        randomSeed(1);
        if (toDrawPhi) {
            shapeDrawer.drawPhi();
        }
        if (toDrawVectors) {
            shapeDrawer.drawVectors();
        }
        shapeDrawer.drawFlyingPoints();
        coordinateSystem.draw();
        shapeDrawer.draw();
        popMatrix();


        // Draw some text (not panned or zoomed!)
        fill(0);
        text("a: zoom in\n" +
                "z: zoom out\n" +
                "p: compute phi(or psi)\n" +
                "l: show phi\n" +
                "v: compute vectors\n" +
                "b: draw vectors\n" +
                "drag mouse to pan", 10, 32);
        shapeDrawer.drawColorManagers();
    }

    // Zoom in and out when the key is pressed
    public void keyPressed() {
        if (key == 'a') {
            zoom += 0.1;
        } else if (key == 'z') {
            zoom -= 0.1;
        } else if (key == 's') {
            solver.solve();
            shapeDrawer.computeColorManagers();
            if (toDrawPhi) {
                shapeDrawer.computePhi();
            }
            if (toDrawVectors) {
                shapeDrawer.computeVectors();
            }
            System.out.println("solved");
        } else if (key == 'p') {
            shapeDrawer.computePhi();
            System.out.println("phi computed");
        } else if (key == 'l') {
            toDrawPhi = !toDrawPhi;
            if (toDrawPhi) {
                shapeDrawer.computePhi();
            }
        } else if (key == 'v') {
            shapeDrawer.computeVectors();
            if (toDrawVectors) {
                shapeDrawer.computeVectors();
            }
            System.out.println("computeVectors computed");
        } else if (key == 'b') {
            toDrawVectors = !toDrawVectors;
        }
        zoom = constrain(zoom, 0, 100);
    }

    // Store the mouse and the previous offset
    public void mousePressed() {
        mouse = new PVector(mouseX, mouseY);
        poffset.set(offset);
    }

    // Calculate the new offset based on change in mouse vs. previous offsey
    public void mouseDragged() {
        offset.x = mouseX - mouse.x + poffset.x;
        offset.y = mouseY - mouse.y + poffset.y;
    }
}
