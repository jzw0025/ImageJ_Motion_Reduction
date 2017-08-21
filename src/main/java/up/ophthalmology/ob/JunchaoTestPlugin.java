/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package up.ophthalmology.ob;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math.ArgumentOutsideDomainException;
import org.apache.commons.math.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;
import org.python.modules.math;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;
import net.imagej.ops.Ops.Copy.Img;
import ij.plugin.filter.GaussianBlur;

/**
 * A template for processing each pixel of either
 * GRAY8, GRAY16, GRAY32 or COLOR_RGB images.
 *
 * @author Junchao Wei 08/15/2017
 */
public class JunchaoTestPlugin implements PlugInFilter {
	protected ImagePlus image;
	private ImageProcessor copy_image; //  ImageProcess 'ip' to create ImagePlus image
	private ImagePlus marking_image;   // ImageProcess 'copy_image' to create ImagePlus marking_image
	// image property members
	private static int width;
	private static int height;
	
	// plugin parameters
	public double value;
	public String name;
	public Map<Integer, Integer> displacement = new HashMap<Integer, Integer>();

	@Override
	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		image = imp;
		return DOES_8G | DOES_16 | DOES_32 | DOES_RGB;
	}
	
	private float mean(float[] input, int start, int center, int end, boolean weighted, double sigma) {
		/* this method calculates the mean of given array*/
		int len0 = end - start + 1;
		double[] weights = new double[len0];
		double var = sigma;
		if (weighted) {
			if (center-start == end-center) {
				for(int i=0; i<len0; i++) {// Centered condition
					weights[i] = 1/Math.sqrt(2*Math.PI)*math.exp(-1.0* math.pow((i+start-center),2)/2*var);
					//System.out.println(weights[i]);
					// Centered condition
				} 
			} else {
				for(int i=0; i<len0; i++) {// Asymmetric condition	
					if (center-start > end-center) {
						if(i+start-center < center - end) { // short radius = center - end
							System.out.println(center);
							weights[i] = 2/Math.sqrt(2*Math.PI)*math.exp(-1.0* math.pow((i+start-center),2)/2*var);
						} else {
							weights[i] = 1/Math.sqrt(2*Math.PI)*math.exp(-1.0* math.pow((i+start-center),2)/2*var);
						}
					}else if(center-start < end-center) {
						if(i+start-center > center - start) { // short radius = center - end
							System.out.println(center);
							weights[i] = 2/Math.sqrt(2*Math.PI)*math.exp(-1.0* math.pow((i+start-center),2)/2*var);
						} else {
							weights[i] = 1/Math.sqrt(2*Math.PI)*math.exp(-1.0* math.pow((i+start-center),2)/2*var);
						}
					}
					else {
						weights[i] = 1/Math.sqrt(2*Math.PI)*math.exp(-1.0* math.pow((i+start-center),2)/2*var);
					}
				}
				
			}		

		}else {
			for(int i=0; i<len0; i++) {
				weights[i] = 1.0/len0; // if not weighted, the weighted is 1/N
			}
		}

		float sum = 0;
		for(int i=start; i<end+1; i++) {
			sum += weights[i-start]*input[i]; //convolve the weights with input;
		}
		return sum;
	}
	
	public double[] linearInterp(double[] x, double[] y, double[] xi) throws ArgumentOutsideDomainException {
		   // return linear interpolation of (x,y) on xi
		   LinearInterpolator li = new LinearInterpolator();
		   PolynomialSplineFunction psf = li.interpolate(x, y);
		   double[] yi = new double[xi.length];
		   for(int i=0; i<xi.length;i++) {
			   yi[i] = psf.value(xi[i]);
		   }
		   return yi;
		}
	
	private float[] double2float(double[] input) {
		float[] array = new float[input.length];
		for(int i=0; i<input.length; i++) {
			array[i] = (float) input[i];
		}
		System.out.println("converted input double array to float array...");
		return array;
	}
	
	private double[] float2double(float[] input) {
		double[] array = new double[input.length];
		for(int i=0; i<input.length; i++) {
			array[i] = input[i];
		}
		System.out.println("converted input float array to double array...");
		return array;
	}
	
	@Override
	public void run(ImageProcessor ip) {
		// get width and height
		width = ip.getWidth();
		height = ip.getHeight();
		
		copy_image = ip.duplicate(); // copy the image for marking the motion.
		
		GaussianBlur gb = new GaussianBlur();
		gb.blurGaussian(copy_image, 1, 1, 0.02);
		
        marking_image = new ImagePlus("Marking the Motion Image", copy_image);
        marking_image.show();

		if (showDialog()) {
			System.out.println(ip);
			//process(ip); // this process the single image slice
			//process(image); // this process the image stacks
			new WaitForUserDialog("User Input Required","Select the Polygon Line, then click OK.").show();
			
			Roi roi = marking_image.getRoi(); // get ROI from marking image
			
			if (!(roi!=null && roi.getType()==Roi.POLYLINE))
				{IJ.error("Straight line selection required."); return;}
			
			FloatPolygon fp_non_inter = ((PolygonRoi) roi).getFloatPolygon();
			
			float[] grid_xc = fp_non_inter.xpoints;
			float[] grid_yc = fp_non_inter.ypoints;
			float[] grid_yc_mean = new float[grid_yc.length];
			
			int begin1 = 0;
			int end1 = 0;
			int center1 = 0;
			for(int i=0; i<grid_yc_mean.length; i++) {
				center1 = i;
				begin1 = i - (int) value >= 0 ? i-(int) value:0; // radius selection: left side condition
				end1 = i + (int) value < grid_yc.length ? i + (int) value: grid_yc.length-1; // radius selection: right side condition
				grid_yc_mean[i] = mean(grid_yc, begin1, center1, end1, true, 1.0); // compute the smoothed points based on the input grid points
			}
			double[] grid_yc_mean_disp = new double[grid_yc_mean.length];
			for(int i=0; i<grid_yc_mean.length; i++) {
				grid_yc_mean_disp[i] = grid_yc_mean[i] - grid_yc[i]; // compute the displacement based the grid point y coordinate difference
			}
			
			double[] grid_xc_double = float2double(grid_xc);
			//double[] grid_yc_mean_double = float2double(grid_yc_mean);
			
			double[] pixel_x = new double[width];
			for(int i=0; i<width; i++) { 
				pixel_x[i] = i; // creating linear interpolation bases;
			}
			
			// establishing interpolation...
			double[] pixel_y = new double[pixel_x.length];
 			try {
				pixel_y = linearInterp(grid_xc_double, grid_yc_mean_disp, pixel_x);
			} catch (ArgumentOutsideDomainException e) {
				System.out.println("the input x is out of range..." );
				e.printStackTrace();
			}
 			
 			float [] pixel_x_float = double2float(pixel_x);
 			float [] pixel_y_float = double2float(pixel_y);
			
 			PolygonRoi drawline = new PolygonRoi(grid_xc, grid_yc_mean, Roi.POLYLINE);
			drawline.drawPixels(copy_image); // ploting the smoothed line onto current image.
			marking_image.updateAndDraw();
			
			for(int i=0; i<pixel_x_float.length; i++) {
				int location = Math.round(pixel_x_float[i]);
				int move = Math.round(pixel_y_float[i]);
				displacement.put(location, move); // put pixel displacement into the hashmap
			}
			
			process(image);
			image.updateAndDraw();
		}
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Process pixels");

		// default value is 0.00, 2 digits right of the decimal point
		gd.addNumericField("Smooth Radius", 0.00,  5);
		gd.addStringField("Version", "Test.v1.0");

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		value = gd.getNextNumber();
		name = gd.getNextString();

		return true;
	}

	/**
	 * Process an image.
	 * <p>
	 * Please provide this method even if {@link ij.plugin.filter.PlugInFilter} does require it;
	 * the method {@link ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)} can only
	 * handle 2-dimensional data.
	 * </p>
	 * <p>
	 * If your plugin does not change the pixels in-place, make this method return the results and
	 * change the {@link #setup(java.lang.String, ij.ImagePlus)} method to return also the
	 * <i>DOES_NOTHING</i> flag.
	 * </p>
	 *
	 * @param image the image (possible multi-dimensional)
	 */
	
	public void process(ImagePlus image) {
		// slice numbers start with 1 for historical reasons
		//System.out.println("Image Stack Detected!");
		for (int i = 1; i <= image.getStackSize(); i++)
			process(image.getStack().getProcessor(i));
	}

	// Select processing method depending on image type
	public void process(ImageProcessor ip) {
		int type = image.getType();
		if (type == ImagePlus.GRAY8)
			process( (byte[]) ip.getPixels() );
		else if (type == ImagePlus.GRAY16)
			process( (short[]) ip.getPixels() );
		else if (type == ImagePlus.GRAY32)
			process( (float[]) ip.getPixels() );
		else if (type == ImagePlus.COLOR_RGB)
			process( (int[]) ip.getPixels() );
		else {
			throw new RuntimeException("not supported");
		}
	}

	// processing of GRAY8 images
	public void process(byte[] pixels) {
		byte[] pixels_copy = new byte[width*height];
		for(int i=0;i<width*height; i++) {
			pixels_copy[i] = 0;
		}
		for (int y=10; y < height-10; y++) {
			for (int x=10; x < width-10; x++) {
				// process each pixel of the line
				// example: add 'number' to each pixel
				int pixel_displacement = displacement.get(x);
				pixels_copy[x + y * width] = pixels[x + y * width - pixel_displacement*width];
			}
		}
		for(int k=0; k<width*height; k++) {
			pixels[k] = pixels_copy[k];
		}
	}

	// processing of GRAY16 images
	public void process(short[] pixels) {
		short[] pixels_copy = new short[width*height];
		for(int i=0;i<width*height; i++) {
			pixels_copy[i] = 0;
		}
		for (int y=10; y < height-10; y++) {
			for (int x=10; x < width-10; x++) {
				// process each pixel of the line
				// example: add 'number' to each pixel
				int pixel_displacement = displacement.get(x);
				pixels_copy[x + y * width] = pixels[x + y * width - pixel_displacement*width];
			}
		}
		for(int k=0; k<width*height; k++) {
			pixels[k] = pixels_copy[k];
		}
	}

	// processing of GRAY32 images
	public void process(float[] pixels) {
		float[] pixels_copy = new float[width*height];
		for(int i=0;i<width*height; i++) {
			pixels_copy[i] = 0;
		}
		for (int y=10; y < height-10; y++) {
			for (int x=10; x < width-10; x++) {
				// process each pixel of the line
				// example: add 'number' to each pixel
				int pixel_displacement = displacement.get(x);
				pixels_copy[x + y * width] = pixels[x + y * width - pixel_displacement*width];
			}
		}
		for(int k=0; k<width*height; k++) {
			pixels[k] = pixels_copy[k];
		}
	}

	// processing of COLOR_RGB images
	public void process(int[] pixels) {
		int[] pixels_copy = new int[width*height];
		for(int i=0;i<width*height; i++) {
			pixels_copy[i] = 0;
		}
		for (int y=10; y < height-10; y++) {
			for (int x=10; x < width-10; x++) {
				// process each pixel of the line
				// example: add 'number' to each pixel
				int pixel_displacement = displacement.get(x);
				pixels_copy[x + y * width] = pixels[x + y * width - pixel_displacement*width];
			}
		}
		for(int k=0; k<width*height; k++) {
			pixels[k] = pixels_copy[k];
		}
	}

	public void showAbout() {
		IJ.showMessage("ProcessPixels",
			"a template for processing each pixel of an image"
		);
	}

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads
	 * an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = JunchaoTestPlugin.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring("file:".length(), url.length() - clazz.getName().length() - ".class".length());
		System.setProperty("plugins.dir", pluginsDir);
		
		System.out.println(url);
		// start ImageJ
		new ImageJ();

		// open the Clown sample
		//ImagePlus image = IJ.openImage("http://imagej.net/images/clown.jpg");
		//ImagePlus image = IJ.openImage("/Users/junchaowei/Dropbox/try/acute3_OD_V_3x3_0_0001016_reslice-1.tif");
		
		ImagePlus image = IJ.openImage(); // This class field image is initialized by open a file.
		image.show();
		
		//System.out.println(image);
		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}