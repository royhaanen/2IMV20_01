/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor. 
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    private int rendererMode = 0;
    double stepSize = 3.0;
        
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] >= volume.getDimX() || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
    }


    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }


    void mip(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        double stepSize = Double.parseDouble(panel.stepSize.getText());

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        int[] volumeDimensions = new int[] {volume.getDimX(), volume.getDimY(), volume.getDimZ()};
        int maxDimension = getMax(volumeDimensions);
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                
                int maxVal = 0;
        
                for (double td = 0; td <= 2 * maxDimension; td+=stepSize) { 
                // Optimialisatie wellicht met 
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + td * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + td * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + td * viewVec[2];
                    
                    if (pixelCoord[0] > volume.getDimX() || pixelCoord[1] > volume.getDimY() || pixelCoord[2] > volume.getDimZ()) {
                        break;
                    }

                    int val = getVoxel(pixelCoord);
                    
                    if (val > maxVal) 
                        maxVal = val;
                }
                
                for (int t = 0; t >= - 2 * maxDimension; t--) {
                    
                    double td = 1 * t;
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + td * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + td * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + td * viewVec[2];
                    
                    // Optimization by checking upfront if pixel coordinates are negative and break the loop
                    if (pixelCoord[0] < 0 || pixelCoord[1] < 0 || pixelCoord[2] < 0) {
                        break;
                    }

                    int val = getVoxel(pixelCoord);
                    
                    if (val > maxVal) 
                        maxVal = val;
                }
                
                // Map the intensity to a grey value by linear scaling
//                voxelColor.r = maxVal/max;
//                voxelColor.g = voxelColor.r;
//                voxelColor.b = voxelColor.r;
//                voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
//                // Alternatively, apply the transfer function to obtain a color
                voxelColor = tFunc.getColor(maxVal);

                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        } 
    }
    
    void compositing(double[] viewMatrix, boolean triLinearInterpolation) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        
        double stepSize = Double.parseDouble(panel.stepSize.getText());

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        int[] volumeDimensions = new int[] {volume.getDimX(), volume.getDimY(), volume.getDimZ()};
        int maxDimension = getMax(volumeDimensions);
        
        double threshold = 1.0;
        
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                
                TFColor prevColor = new TFColor();
                TFColor nextColor = new TFColor(); 
                
                for (double t = - 0.5 * maxDimension; t <= 0.5 * maxDimension; t+=stepSize) {
                // Optimization possible by step size
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + t * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + t * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + t * viewVec[2];

                    
                    int val = 0;
                    if (triLinearInterpolation) {
                        val = triLinearInterpolation(pixelCoord);
                    }
                    else {
                        val = getVoxel(pixelCoord);
                    }
                    
                    if ( ( (pixelCoord[0] < volume.getDimX() && pixelCoord[0] >= 0) || (pixelCoord[1] < volume.getDimY() && pixelCoord[1] >= 0) || (pixelCoord[2] < volume.getDimZ() && pixelCoord[2] >= 0) ) && val > threshold) {
                        
                        voxelColor = tFunc.getColor(val);
                        
                        nextColor.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * prevColor.r;
                        nextColor.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * prevColor.g;
                        nextColor.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * prevColor.b;
                        
                        nextColor.a = (1 - voxelColor.a) * prevColor.a;
                        
                        prevColor = nextColor;
                    }
                }
                
                //System.out.println("op: " + (1 - nextColor.a));
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = (1 - nextColor.a) <= 1.0 ? (int) Math.floor((1 - nextColor.a) * 255) : 255;
                int c_red = nextColor.r <= 1.0 ? (int) Math.floor(nextColor.r * 255) : 255;
                int c_green = nextColor.g <= 1.0 ? (int) Math.floor(nextColor.g * 255) : 255;
                int c_blue = nextColor.b <= 1.0 ? (int) Math.floor(nextColor.b * 255) : 255;
                    
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        } 
    }
    
    void twodtransfer(double[] viewMatrix, boolean illumination) {

         // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        double stepSize = Double.parseDouble(panel.stepSize.getText());
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // Illumination variables
        double[] lVec = new double[3]; // L-Vector
        double[] hVec = new double[3]; // H-Vector
        double[] nVec = new double[3]; // N-Vector
        lVec = viewVec; // Light is coming along the view Vector
        hVec = lVec; // Therefore L Vector is equal to H Vector
        double i_a = 0.0; // ambient light
        int alpha = 10;
        double k_spec = 0.2;
        double k_diff = 0.7;
        double k_ambient = 0.7;
        // Assuming White color light
        
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        int[] volumeDimensions = new int[] {volume.getDimX(), volume.getDimY(), volume.getDimZ()};
        int maxDimension = getMax(volumeDimensions);
        
        double threshold = 1.0;
        
        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        
        // Get the GUI user defined variables
        int definedIntensity = tfEditor2D.triangleWidget.baseIntensity;
        double definedRadius = tfEditor2D.triangleWidget.radius;
        TFColor definedColor = tfEditor2D.triangleWidget.color;
        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                
                TFColor prevColor = new TFColor();
                TFColor nextColor = new TFColor(); 
                
                for (double t = - 0.5 * maxDimension; t <= 0.5 * maxDimension; t+=stepSize) {
                // Optimization possible by step size
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0] + t * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1] + t * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2] + t * viewVec[2];

                    // Use the line below for non-interpolated values
                    int val = getVoxel(pixelCoord);
                    
                    // Use the line below for tri-linear interpolated values
                    //int val = triLinearInterpolation(pixelCoord);
                    
                    if ( ( (pixelCoord[0] < volume.getDimX() && pixelCoord[0] >= 0) || (pixelCoord[1] < volume.getDimY() && pixelCoord[1] >= 0) || (pixelCoord[2] < volume.getDimZ() && pixelCoord[2] >= 0) ) && val > threshold) {
                        
                        VoxelGradient voxGrad = gradients.getGradient((int)Math.floor(pixelCoord[0]), (int)Math.floor(pixelCoord[1]), (int)Math.floor(pixelCoord[2])); 
                        
                        if (val == definedIntensity && voxGrad.mag == 0) {
                            voxelColor.a = definedColor.a * 1.0;
                        }
                        else if (voxGrad.mag > 0.0 && ((val - definedRadius * voxGrad.mag) <= definedIntensity) && ((val + definedRadius * voxGrad.mag) >= definedIntensity) )  {
                            voxelColor.a = definedColor.a * (1.0 - (1 / definedRadius) * (Math.abs((definedIntensity - val)/ voxGrad.mag)));
                        }
                        else 
                            voxelColor.a = 0.0;
                        
                        if (illumination) {
                            if (voxGrad.mag > 0.0 && voxelColor.a > 0.0) {
                                // Filling N-Vector:
                                nVec[0] = voxGrad.x / voxGrad.mag;
                                nVec[1] = voxGrad.y / voxGrad.mag;
                                nVec[2] = voxGrad.z / voxGrad.mag;

                                // Computing required dot products
                                double l_dot_n = VectorMath.dotproduct(lVec, nVec);
                                double n_dot_h = VectorMath.dotproduct(nVec, hVec);

                                if (l_dot_n > 0 && n_dot_h > 0 ) {
                                    voxelColor.r = i_a + (definedColor.r * k_diff * l_dot_n) + k_spec * Math.pow(n_dot_h, alpha);
                                    voxelColor.g = i_a + (definedColor.g * k_diff * l_dot_n) + k_spec * Math.pow(n_dot_h, alpha);
                                    voxelColor.b = i_a + (definedColor.b * k_diff * l_dot_n) + k_spec * Math.pow(n_dot_h, alpha);
                                }
                            }
                            
                            nextColor.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * prevColor.r;
                            nextColor.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * prevColor.g;
                            nextColor.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * prevColor.b;
                        }
                         
                        else {
                            nextColor.r = voxelColor.a * definedColor.r + (1 - voxelColor.a) * prevColor.r;
                            nextColor.g = voxelColor.a * definedColor.g + (1 - voxelColor.a) * prevColor.g;
                            nextColor.b = voxelColor.a * definedColor.b + (1 - voxelColor.a) * prevColor.b;
                        }
                        
                        nextColor.a = (1 - voxelColor.a) * prevColor.a;
                        
                        prevColor = nextColor;
                        
                    }
                    
                }
                        
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = (1 - nextColor.a) <= 1.0 ? (int) Math.floor((1 - nextColor.a) * 255) : 255;
                int c_red = nextColor.r <= 1.0 ? (int) Math.floor(nextColor.r * 255) : 255;
                int c_green = nextColor.g <= 1.0 ? (int) Math.floor(nextColor.g * 255) : 255;
                int c_blue = nextColor.b <= 1.0 ? (int) Math.floor(nextColor.b * 255) : 255;
                    
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        } 

    }
    
    public static double linearInt(double x, double x0, double x1, double v0, double v1) {
        
        double v = ((x1 - x) / (x1 - x0)) * v0 + ((x - x0) / (x1 - x0)) * v1;
        
        return v;
    }

    public static short triLinearInt(double x, double y, double z, double q000, double q001, double q010, double q011, double q100, double q101, double q110, double q111, double x1, double x2, double y1, double y2, double z1, double z2) {
      // x =  pixelCoord[0], y = pixelCoord[1], z = pixelCoord[2] 
      // q000 = value of (floor(x), floor(y), floor(z))...
        
      double x00 = linearInt(x, x1, x2, q000, q100);
      double x10 = linearInt(x, x1, x2, q010, q110);
      double x01 = linearInt(x, x1, x2, q001, q101);
      double x11 = linearInt(x, x1, x2, q011, q111);
      double r0 = linearInt(y, y1, y2, x00, x01);
      double r1 = linearInt(y, y1, y2, x10, x11);

      return (short)Math.round(linearInt(z, z1, z2, r0, r1));
    }
    
    public short triLinearInterpolation (double[] pixelCoord) {
        
        short val = 0;
        
        double x = pixelCoord[0];
        double y = pixelCoord[1];
        double z = pixelCoord[2];
        
        double q000[] = {Math.floor(x), Math.floor(y), Math.floor(z)};
        double q100[] = {Math.floor(x)+1.0, Math.floor(y), Math.floor(z)};
        double q010[] = {Math.floor(x), Math.floor(y)+1.0, Math.floor(z)};
        double q110[] = {Math.floor(x)+1.0, Math.floor(y)+1.0, Math.floor(z)};
        double q001[] = {Math.floor(x), Math.floor(y), Math.floor(z)+1.0};
        double q101[] = {Math.floor(x)+1.0, Math.floor(y), Math.floor(z)+1.0};
        double q011[] = {Math.floor(x), Math.floor(y)+1.0, Math.floor(z)+1.0};
        double q111[] = {Math.floor(x)+1.0, Math.floor(y)+1.0, Math.floor(z)+1.0};
        
        val = triLinearInt(x, y, z, (double)getVoxel(q000), (double)getVoxel(q001), (double)getVoxel(q010), (double)getVoxel(q011), (double)getVoxel(q100), (double)getVoxel(q101), (double)getVoxel(q110), (double)getVoxel(q111), Math.floor(x), Math.floor(x)+1.0, Math.floor(y), Math.floor(y)+1.0, Math.floor(z), Math.floor(z)+1.0);
        
        return val;
    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }
    
    // Method for getting the maximum value
    public static int getMax(int[] inputArray){ 
    
        int maxValue = inputArray[0]; 
        
        for(int i=1;i < inputArray.length;i++){ 
         
            if(inputArray[i] > maxValue){ 
            
                maxValue = inputArray[i]; 
            } 
        } 
    
        return maxValue; 
  }

    @Override
    public void visualize(GL2 gl) {

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);
        
        long startTime = System.currentTimeMillis();
        
        // rendererModes: 0 = slicer, 1 = mip, 2 = compositing, 3 = 2dtransfer
        switch (rendererMode) {
            case 0:     panel.shadingCheckbox.setSelected(false);
                        slicer(viewMatrix);
                        break;
            case 1:     panel.shadingCheckbox.setSelected(false);
                        mip(viewMatrix);
                        break;
            case 2:     panel.shadingCheckbox.setSelected(false);
                        panel.compositingButton.setSelected(true);
                        compositing(viewMatrix, panel.triLinearCheckbox.isSelected());
                        break;
            case 3:     twodtransfer(viewMatrix, panel.shadingCheckbox.isSelected());
                        break;
            case 4:     panel.tf2dButton.setSelected(true);
                        twodtransfer(viewMatrix, panel.shadingCheckbox.isSelected());
                        break;
            default:    System.out.println("Invalid renderer method.");
                        break;
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
    
    public void setMode(int mode) {
        rendererMode = mode;
    }
}
