package cseis.seaview;

import java.util.ArrayList;
import cseis.seisdisp.csWellPathAttr;
import cseis.seisdisp.csWellPathOverlayAttr;

public interface csIWellPathDialogListener {
  public void updateWellPaths( ArrayList<csWellPathAttr> wellPathAttrList, csWellPathOverlayAttr wellPathOverlayAttr );
  public void updateWellPathDirectory( String directory );
}
