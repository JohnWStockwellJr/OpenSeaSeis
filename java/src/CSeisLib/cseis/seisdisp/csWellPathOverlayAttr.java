package cseis.seisdisp;

public class csWellPathOverlayAttr {
  public static final int MIN_TRANSPARENCY = 50;
  public static final int MAX_TRANSPARENCY = 255;
  public static final int MIN_SHADOW_DEPTH = 1;
  public static final int MAX_SHADOW_DEPTH = 250;
  /// Transparency of well path plot (50-255)
  public int     transparency;
  /// Physical distance in [m] from seismic section to which well path is plotted
  public int     maxDisplayDistance;
  /// Size of line plot
  public int     lineSize;
  /// Show well path shadow
  public boolean showShadow;
  /// A number between 10 (small depth) to 1000 (steep depth). 0: 'Sun' is shining from front. infinity: 'Sun' is shining from top.
  public int shadowDepth;
  /// Trace header name in seismic containing X coordinate
  public String hdrName_binx;
  /// Trace header name in seismic containing Y coordinate
  public String hdrName_biny;
  //  public double shadowDirectionAngle;
  
  public csWellPathOverlayAttr() {
    transparency = 255;
    maxDisplayDistance = 5000;
    showShadow  = false;
    shadowDepth = 100;
    lineSize    = 4;
    hdrName_binx = "bin_x";
    hdrName_biny = "bin_y";
  }
  public csWellPathOverlayAttr( csWellPathOverlayAttr attr ) {
    transparency = attr.transparency;
    maxDisplayDistance = attr.maxDisplayDistance;
    showShadow  = attr.showShadow;
    shadowDepth = attr.shadowDepth;
    lineSize    = attr.lineSize;
    hdrName_binx = attr.hdrName_binx;
    hdrName_biny = attr.hdrName_biny;
  }
}
