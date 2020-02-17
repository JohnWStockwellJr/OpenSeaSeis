/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csColorBarPanel;
import cseis.general.csColorMap;
import cseis.seisdisp.csAnnotationAttributes;
import cseis.seisdisp.csSampleInfo;
import cseis.seisdisp.csSeisDispSettings;
import cseis.seisdisp.csSeisPane;
import cseis.seisdisp.csSeisView;
import java.awt.BorderLayout;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import javax.swing.BorderFactory;
import javax.swing.JScrollBar;
import javax.swing.JSplitPane;
import cseis.seisdisp.csWellPathAttr;
import cseis.seisdisp.csWellPathOverlayAttr;

/**
 * The csSeisPaneManager manages all seismic bundles contained in SeaView.<br>
 * Helps relay method calls between encompassing classes (SeaView) and seismic bundles.
 * @author 2013 Felipe Punto
 */
// TODO: Most methods are really superfluous since they simply call similar method in seismic bundle.
// Potentially eliminate this class and move its little genuine functionality to csSeisPaneBundle or SeaView.
public class csSeisPaneManager {
  private ArrayList<csSeisPaneBundle> myBundleList;
  private SeaView mySeaView;
  private boolean myIsSyncUpdating;
  private boolean myIsSyncDispSettings; // Sync display settings across sync'ed seimsic bundles
  private int myCrosshairStatus;
  private static final int CROSSHAIR_OFF = 1;
  private static final int CROSSHAIR_ON_NOSNAP = 2;
  private static final int CROSSHAIR_ON_SNAP  = 3;

  public csSeisPaneManager( SeaView seaview, csAnnotationAttributes annAttr ) {
    myBundleList = new ArrayList<csSeisPaneBundle>();
    mySeaView = seaview;
    myIsSyncUpdating = false;
    myCrosshairStatus = csSeisPaneManager.CROSSHAIR_OFF;
    myIsSyncDispSettings = false;
  }
  public void updateSyncDispSettings( boolean doSync ) {
    myIsSyncDispSettings = doSync;
  }
  public csSeisPaneBundle add( csSeisView seisViewIn, csTraceSelectionParam traceSelectionParam, cseis.general.csFilename filename, int fileFormat ) {
    csSeisPaneBundle bundle = new csSeisPaneBundle( mySeaView, this, seisViewIn, traceSelectionParam, filename, fileFormat );

    bundle.seisPane = new csSeisPane( bundle.seisView, mySeaView.getAnnotationAttributes() );
    bundle.seisPane.setCrosshair( myCrosshairStatus != csSeisPaneManager.CROSSHAIR_OFF, myCrosshairStatus == csSeisPaneManager.CROSSHAIR_ON_SNAP );

    bundle.splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, false, null, bundle.seisPane );
    bundle.splitPane.setDividerSize( 0 );
    bundle.splitPane.setResizeWeight( 0 );
    bundle.splitPane.setBorder( BorderFactory.createEmptyBorder() );
    bundle.add( bundle.splitPane,BorderLayout.CENTER );
    bundle.popupMenu = new csSeaViewPopupMenu( bundle );
    bundle.seisView.setPopupMenu( bundle.popupMenu );
    bundle.graphPanel = null;
    bundle.colorBarPanel = null;
    bundle.addSeisPaneBundleMouseModeListener( bundle.popupMenu );
    myBundleList.add(bundle);

    bundle.addSeisPaneBundleScrollListener( new csISeisPaneBundleScrollListener() {
      @Override
      public void horzScrollChanged( csSeisPaneBundle bundle, int value ) {
        if( bundle.isSynced() ) {
          if( myIsSyncUpdating ) return;
          myIsSyncUpdating = true;
          for( int i = 0; i < getNumBundles(); i++ ) {
            csSeisPaneBundle bundle2 = myBundleList.get(i);
            if( !bundle.equals(bundle2) && bundle2.isSynced() ) {
              JScrollBar bar = bundle2.seisPane.getHorizontalScrollBar();
              if( value <= bar.getMaximum() ) bar.setValue(value);
            }
          }
          myIsSyncUpdating = false;
        }
      }
      @Override
      public void vertScrollChanged( csSeisPaneBundle bundle, int value ) {
        if( bundle.isSynced() ) {
          if( myIsSyncUpdating ) return;
          myIsSyncUpdating = true;
          for( int i = 0; i < getNumBundles(); i++ ) {
            csSeisPaneBundle bundle2 = myBundleList.get(i);
            if( !bundle.equals(bundle2) && bundle2.isSynced() ) {
              JScrollBar bar = bundle2.seisPane.getVerticalScrollBar();
              if( value <= bar.getMaximum() ) bar.setValue(value);
            }
          }
          myIsSyncUpdating = false;
        }
      }
    });
    bundle.addSeisPaneBundleListener( new csISeisPaneBundleListener() {
      @Override
      public void updateBundleDisplaySettings( csSeisPaneBundle bundle ) {
        if( !myIsSyncDispSettings ) return;
        if( bundle.isSynced() ) {
          if( myIsSyncUpdating ) return;
          myIsSyncUpdating = true;
          csSeisDispSettings ds = bundle.seisView.getDispSettings();
          for( int i = 0; i < getNumBundles(); i++ ) {
            csSeisPaneBundle bundle2 = myBundleList.get(i);
            if( !bundle.equals(bundle2) && bundle2.isSynced() ) {
              bundle2.seisView.updateDispSettings(ds);
            }
          }
          myIsSyncUpdating = false;
        }
      }
      @Override
      public void updateBundleDisplayScalar( csSeisPaneBundle bundle ) {
        // Nothing to be done
      }
      @Override
      public void updateBundleSampleInfo( csSeisPaneBundle bundle, csSeisBundleSampleInfo sampleInfo ) {
                  // Nothing
      }
    });
   
    return bundle;
  }
  public csAnnotationAttributes getAnnotationAttributes() {
    return mySeaView.getAnnotationAttributes();
  }
  public void removeBundle( csSeisPaneBundle bundle ) {
    bundle.close();
    myBundleList.remove( bundle );
  }
  public int getNumBundles() {
    return myBundleList.size();
  }
  public csSeisPaneBundle getBundle( int index ) {
    return myBundleList.get(index);
  }
  public void showColorBar( boolean doShow ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      csSeisPaneBundle bundle = myBundleList.get(i);
      showColorBar( doShow, bundle );
    } // END for bundleList
  }
  public void showColorBar( boolean doShow, csSeisPaneBundle bundle ) {
    if( bundle.colorBarPanel == null ) {
      bundle.colorBarPanel = new csColorBarPanel( new csColorMap() );
      bundle.setColorBar( bundle.seisView.getDispSettings() );
    }
    if( doShow ) {
      bundle.add( bundle.colorBarPanel, BorderLayout.EAST );
    }
    else {
      bundle.remove( bundle.colorBarPanel );
      bundle.colorBarPanel = null;
    }
  }
  public void setShowFilenameInView( boolean showFilename, csSeisPaneBundle activeBundle ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      csSeisPaneBundle bundle = myBundleList.get(i);
      bundle.setShowFilenameInView( showFilename );
      if( bundle.equals(activeBundle) ) bundle.seisView.repaint();
    } // END for bundleList
  }
  public void updateWellPaths( ArrayList<csWellPathAttr> wellPathAttrList, csWellPathOverlayAttr wellPathOverlayAttr ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      csSeisPaneBundle bundle = myBundleList.get(i);
      bundle.updateWellPaths( wellPathAttrList, wellPathOverlayAttr );
    } // END for bundleList
  }
  //-------------------------------------------------------------
  public void zoom( int zoomType, int zoomMode, csSeisPaneBundle bundle ) {
    if( bundle != null ) {
      bundle.seisPane.zoom( zoomType, zoomMode );
    }
  }
  public void setColorBits( int colorMapType ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      myBundleList.get(i).seisView.setColorBits(colorMapType);
    }
  }
  public void highlightTrace( boolean doHighlightTrace ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      highlightTrace( doHighlightTrace, myBundleList.get(i) );
    }
  }
  public void highlightTrace( boolean doHighlightTrace, csSeisPaneBundle bundle ) {
    bundle.seisPane.setAutoHighlightTrace( doHighlightTrace );
  }
  public void updateScalar( float scalar, csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.updateScalar( scalar );
  }
  public void incScaling( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.incScaling();
  }
  public void decScaling( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.decScaling();
  }
  public void showSettingsDialog( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.showSettingsDialog();
  }
  public void showTraceMonitor( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.showTraceMonitor();
  }
  public void showHeaderMonitor( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.showHeaderMonitor();
  }
  public void showAnnotationDialog( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.showAnnotationDialog();
  }
  public void showOverlayDialog( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.showOverlayDialog();
  }
//  public void setPaintMode( csSeisPaneBundle bundle, boolean doSet ) {
//    if( bundle != null ) bundle.setPaintMode( doSet );
//  }
//  public void setPickMode( csSeisPaneBundle bundle, boolean doSet ) {
//    if( bundle != null ) bundle.setPickMode( doSet );
//  }
  //---------------------------------------------------------------
  // Graph panel
  //
  public void resetGraphSize() {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      myBundleList.get(i).resetGraphSize();
    }
  }
//  public void setGraphPanel( boolean doShowGraph, csSeisPaneBundle bundle ) {
 //   if( bundle != null ) bundle.setGraphPanel(doShowGraph);
//  }
  public void updateGraph( csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.updateGraph();
  }
  public void showGraphPanel( boolean doShow, csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.showGraph(doShow);
  }
//---------------------------------------------------------------
  public void updateDisplaySettings( csSeisDispSettings settings, csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.seisView.updateDispSettings(settings);
  }
  public void retrieveAnnotationSettings( csAnnotationAttributes annAttr, csSeisPaneBundle bundle ) {
    if( bundle != null ) bundle.retrieveAnnotationSettings( annAttr );
  }
  public void setColorBar( csSeisDispSettings settings ) {
  }
  public void setLogScaleY( boolean setLogScaleY, csSeisPaneBundle bundle ) {
    if( bundle != null ) {
      csSeisDispSettings settings = bundle.seisView.getDispSettings();
      settings.isLogScale = setLogScaleY;
      bundle.seisView.updateDispSettings( settings );
      bundle.seisPane.zoom( settings.zoomVert, settings.zoomHorz );
    }
  }
  public void setMouseMode( int mouseMode ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      myBundleList.get(i).setMouseMode( mouseMode );
    }
  }
  public void setCrosshair( boolean setOn, boolean setSnap ) {
    for( int i = 0; i < myBundleList.size(); i++ ) {
      myBundleList.get(i).seisPane.setCrosshair( setOn, setSnap );
    }
    if( setOn ) {
      myCrosshairStatus = setSnap ? csSeisPaneManager.CROSSHAIR_ON_SNAP : csSeisPaneManager.CROSSHAIR_ON_NOSNAP;
    }
    else {
      myCrosshairStatus = csSeisPaneManager.CROSSHAIR_OFF;
    }
  }
  public boolean isCrosshairSet() {
    return( myCrosshairStatus != csSeisPaneManager.CROSSHAIR_OFF );
  }
  public void setCrosshairPosition( csSeisPaneBundle bundle, csSampleInfo sampleInfo ) {
    for( int i = 0; i < getNumBundles(); i++ ) {
      csSeisPaneBundle bundle2 = getBundle(i);
      if( !bundle.equals(bundle2) && bundle2.isSynced() ) {
        bundle2.seisView.setCrosshairPosition( sampleInfo );
      }
    }
  }
  public void mouseExitedBundle( csSeisPaneBundle bundle ) {
    if( isCrosshairSet() ) {
      for( int i = 0; i < getNumBundles(); i++ ) {
        csSeisPaneBundle bundle2 = getBundle(i);
        if( !bundle.equals(bundle2) && bundle2.isSynced() ) {
         bundle2.seisView.removeCrosshair();
        }
      }
    }
  }
  public BufferedImage createSnapShot( boolean includeSideLabels, csSeisPaneBundle bundle ) {
    if( bundle == null ) return null;
    Rectangle rect = null;
    if( includeSideLabels ) {
      rect = bundle.seisPane.getVisibleRect();
    }
    else {
      rect = bundle.seisView.getVisibleRect();
    }
    BufferedImage image = new BufferedImage(rect.width, rect.height, BufferedImage.TYPE_INT_RGB);
    if( includeSideLabels ) {
      bundle.seisPane.paintAll(image.createGraphics());
    }
    else {
      bundle.seisView.paintAll(image.createGraphics());
    }
    return image;
  }
}

