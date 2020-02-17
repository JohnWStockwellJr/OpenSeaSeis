/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.velocity;

import cseis.seis.csHeader;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seisdisp.csHorizonAttr;
import cseis.seisdisp.csIKeyListener;
import cseis.seisdisp.csISeisOverlay;
import cseis.seisdisp.csSampleInfo;
import cseis.seisdisp.csSeisView;
import java.awt.BasicStroke;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeMap;
import javax.swing.SwingUtilities;

/**
 * Velocity pick overlay over seismic semblance plot.
 * - Display velocity-time functions
 * - Velocity picking functionality
 * 
 */
public class csVelPickOverlay implements csISeisOverlay, MouseMotionListener, MouseListener, csIKeyListener {
  private final csSeisView mySeisView;

  /// List of velocity function maps, one for each velocity field. A velocity function map consists of pairs of ensemble keys and a time-velocity function
  private final ArrayList< TreeMap<csEnsemble,csVelFunction> >  myVelFieldList;
  /// Display attributes for all velocity fields
  private final ArrayList<csHorizonAttr> myDisplayAttrList;
  /// Information on ensemble headers incl velocity trace header
  private final csVelEnsembleInfo myEnsInfo;

  /// Bundle of currently displayed ensembles and actively picked velocity functions
  private final ArrayList<csVelPickBundle> myCurrentActiveBundleList;

  /// true: Active velocity picking is ongoing
  private boolean myIsActive;
  /// Index of velocity field which is currently active
  private int myActiveVelFieldIndex = 0;
  private boolean myIsBundleLock;

  /// Data set sample interva;
  private final float mySampleInt;
  /// Data set number of samples
  private int myNumSamples;
  /// Data set number of currently displayed traces
  private int myNumTraces;
  
  public csVelPickOverlay( csSeisView seisView, csVelEnsembleInfo ensInfo ) {
    mySeisView  = seisView;
    myEnsInfo = ensInfo;
    myIsBundleLock = false;
    myVelFieldList = new ArrayList< TreeMap<csEnsemble,csVelFunction> >(1);
    myDisplayAttrList = new ArrayList<csHorizonAttr>(1);
    mySampleInt = (float)seisView.getSampleInt();
    
    myIsActive = false;
    myActiveVelFieldIndex = -1;
    
    myCurrentActiveBundleList = new ArrayList<csVelPickBundle>();
  }
  public void activate() {
    if( myIsActive ) return;
    mySeisView.addMouseListener( this );
    mySeisView.addMouseMotionListener( this );
    mySeisView.getEventHandler().addKeyListener( this );
    myIsActive = true;
    copyActiveVelFieldToBundle();
  }
  public void deactivate() {
    if( !myIsActive ) return;
    copyActiveVelFieldFromBundle( true );
    mySeisView.removeMouseListener( this );
    mySeisView.removeMouseMotionListener( this );
    mySeisView.getEventHandler().removeKeyListener( this );
    myIsActive = false;
  }
  public int getNumVelFields() {
    return myVelFieldList.size();
  }
  public csVelEnsembleInfo getEnsembleInfo() {
    return myEnsInfo;
  }
  public TreeMap<csEnsemble,csVelFunction> getVelField( int velFieldID ) {
    int index = retrieveVelFieldIndex( velFieldID );
    if( index >= 0 ) {
      if( myIsActive ) copyActiveVelFieldFromBundle( false );
      return myVelFieldList.get(index);
    }
    return null;
  }
  public void setActiveObject( int velFieldID ) {
    int index = retrieveVelFieldIndex( velFieldID );
    if( index >= 0 ) {
      if( myActiveVelFieldIndex == index ) return;
      if( myIsActive && myActiveVelFieldIndex >= 0 ) copyActiveVelFieldFromBundle(true);
      myActiveVelFieldIndex = index;
      copyActiveVelFieldToBundle();
    }
  }
  public void removeVelField( int velFieldID ) {
    int index = retrieveVelFieldIndex( velFieldID );
    if( index >= 0 ) {
      if( myIsActive && myActiveVelFieldIndex == index ) {
        // Clear bundle velFunc...: (copying is really not needed)
        copyActiveVelFieldFromBundle( true );
      }
      myActiveVelFieldIndex = -1;
      myDisplayAttrList.remove(index);
      myVelFieldList.remove(index);
    }
  }
  public void addVelField( csHorizonAttr attr_in ) {
    addVelField( attr_in,  new java.util.ArrayList<csVelFunction>() );
  }
  public void addVelField( csHorizonAttr attr_in, java.util.List<csVelFunction> velFuncList ) {
    TreeMap<csEnsemble,csVelFunction> map = new TreeMap<csEnsemble,csVelFunction>( new Comparator<Object>() {
      @Override
      public int compare( Object o1, Object o2 ) {
        return ((csEnsemble)o1).compareTo( (csEnsemble)o2 );
      }
    });
    for( csVelFunction velFunc : velFuncList ) {
      map.put( velFunc.getEns(), velFunc );
    }
    myDisplayAttrList.add(attr_in);
    myVelFieldList.add( map );
  }
  private int retrieveVelFieldIndex( int velFieldID ) {
    for( int i = 0; i < myDisplayAttrList.size(); i++ ) {
      if( myDisplayAttrList.get(i).getID() == velFieldID ) {
        return i;
      }
    }
    return -1;
  }
  public void refresh() {
    mySeisView.repaint();
  }
  public void updateDisplayAttr( csHorizonAttr attr_in ) {
    int index = retrieveVelFieldIndex( attr_in.getID() );
    if( index >= 0 ) {
      myDisplayAttrList.set( index, attr_in );
      refresh();
    }
  }
  public void traceBufferChanged( csISeismicTraceBuffer traceBuffer ) throws Exception {
    if( myIsActive ) copyActiveVelFieldFromBundle( false );
    
    myNumSamples = traceBuffer.numSamples();
    myNumTraces = traceBuffer.numTraces();
    int counter = 0;
    while( myIsBundleLock ) {
      Thread.sleep(1);
      counter += 1;
      if( counter == 1000 ) {
        break;
//        throw( new Exception("csVelPickOverlay:traceBufferChanged is unexpectedly waiting for draw function to complete... System overload?") );
      }
    }
    // Next, create bundle list with all currently displayed ensembles. Velocity functions in bundle are loaded later
    myIsBundleLock = true;
    myCurrentActiveBundleList.clear();

    // (1) Compute velocity increment in input data
    if( myNumTraces < 2 ) {
      return; // ???
    }
    float vel1 = traceBuffer.headerValues(0)[myEnsInfo.hdrIndexVel].floatValue();
    float vel2 = traceBuffer.headerValues(1)[myEnsInfo.hdrIndexVel].floatValue();
    if( vel2 <= vel1 ) {
      if( myNumTraces < 3 ) {
        myIsBundleLock = false;
        return; // ???
      }
      float vel3 = traceBuffer.headerValues(2)[myEnsInfo.hdrIndexVel].floatValue();
      if( vel3 <= vel2 ) {
        myIsBundleLock = false;
        throw( new Exception("Selected velocity trace header (" + myEnsInfo.hdrNameVel + ") does not appear to contain valid velocities.\n" +
          "Traces 1, 2 and 3 of current gather contain velocity header values of " + vel1 + ", " + vel2 + " and " + vel3 ) );
      }
      vel1 = vel2;
      vel2 = vel3;
    }
    float velIncPerTrace = vel2 - vel1;
    if( myEnsInfo.getDim() == csEnsemble.DIM_1D ) {
      csHeader[] hdrValues = traceBuffer.headerValues(0);
      float vel = hdrValues[myEnsInfo.hdrIndexVel].floatValue();
      myCurrentActiveBundleList.add( new csVelPickBundle( new csEnsemble(), 0, vel, velIncPerTrace ) );
    }
    else {
      for( int itrc = 0; itrc < myNumTraces; itrc++ ) {
        csHeader[] hdrValues = traceBuffer.headerValues(itrc);
        float vel = hdrValues[myEnsInfo.hdrIndexVel].floatValue();
        csEnsemble ens;
        if( myEnsInfo.getDim() == csEnsemble.DIM_2D ) {
          ens = new csEnsemble( hdrValues[myEnsInfo.hdrIndexEns1].intValue() );
        }
        else {
          ens = new csEnsemble( hdrValues[myEnsInfo.hdrIndexEns1].intValue(), hdrValues[myEnsInfo.hdrIndexEns2].intValue() );
        }

        boolean found = false;
        if( !myCurrentActiveBundleList.isEmpty() ) {
          for( csVelPickBundle bundle : myCurrentActiveBundleList ) {
            if( ens.equals( bundle.ens ) ) {
              found = true;
              break;
            }
          }
        }
        if( !found || myCurrentActiveBundleList.isEmpty() ) {
          myCurrentActiveBundleList.add( new csVelPickBundle( ens, itrc, vel, velIncPerTrace ) );
        }
      }
    }
    copyActiveVelFieldToBundle();
    myIsBundleLock = false;
  }
  /**
   * Load (=copy) new active velocity function from main velocity field list into local bundle object, for editing
   */
  private void copyActiveVelFieldToBundle() {
    if( !myIsActive || myActiveVelFieldIndex < 0 ) return;
    TreeMap<csEnsemble,csVelFunction> map = myVelFieldList.get(myActiveVelFieldIndex);
    for( csVelPickBundle bundle : myCurrentActiveBundleList ) {
      csVelFunction velFunc = map.get(bundle.ens);
      if( velFunc != null ) bundle.velFunc = velFunc;
    }
  }
  /**
   * Copy currently active/edited velocity field to main list
   */
  private void copyActiveVelFieldFromBundle( boolean clearVelFunc ) {
    if( myActiveVelFieldIndex < 0 ) return;
    TreeMap<csEnsemble,csVelFunction> map = myVelFieldList.get(myActiveVelFieldIndex);
    for( csVelPickBundle bundle : myCurrentActiveBundleList ) {
      if( bundle.velFunc.numPicks() != 0 ) {
        map.put( bundle.ens, bundle.velFunc );
        if( clearVelFunc ) bundle.velFunc = new csVelFunction( bundle.ens );
      }
      else if( map.containsKey(bundle.ens) ) {
        map.remove( bundle.ens );
      }
    }
  }
  private float velToTrace( float velocity, csEnsemble ens ) {
    if( myCurrentActiveBundleList.isEmpty() ) return -1;
    csVelPickBundle bundlePrev = myCurrentActiveBundleList.get(0);
    if( myCurrentActiveBundleList.size() > 1 ) {
      int ibundle = 1;
      while( !bundlePrev.ens.equals(ens) && ibundle < myCurrentActiveBundleList.size() ) {
        bundlePrev = myCurrentActiveBundleList.get(ibundle++);
      }
    }
    if( !bundlePrev.ens.equals(ens) ) return -1;
    return( bundlePrev.indexOfTrace1 + (velocity - bundlePrev.velOfTrace1) / bundlePrev.velInc );
  }
  private csVelPickBundle traceToBundle( float traceIndex ) {
    if( myCurrentActiveBundleList.isEmpty() ) return null;
    csVelPickBundle bundlePrev = myCurrentActiveBundleList.get(0);
    if( myCurrentActiveBundleList.size() > 1 ) {
      for( int i = 1; i < myCurrentActiveBundleList.size(); i++ ) {
        csVelPickBundle bundle = myCurrentActiveBundleList.get(i);
        if( bundle.indexOfTrace1 > traceIndex ) break;
        bundlePrev = bundle;
      }
    }
    return( bundlePrev );
  }
//  private float traceToVel( float traceIndex ) {
//    if( myCurrentActiveBundleList.isEmpty() ) return 0;
//    csVelPickBundle bundlePrev = myCurrentActiveBundleList.get(0);
//    if( myCurrentActiveBundleList.size() > 1 ) {
//      for( int i = 1; i < myCurrentActiveBundleList.size(); i++ ) {
//        csVelPickBundle bundle = myCurrentActiveBundleList.get(i);
//        if( bundle.indexOfTrace1 > traceIndex ) break;
//        bundlePrev = bundle;
//      }
//    }
//    return( bundlePrev.velOfTrace1 + bundlePrev.velInc * (traceIndex-bundlePrev.indexOfTrace1) );
//  }
//  private csVelFunction traceToVelFunc( float traceIndex ) {
//    if( myCurrentActiveBundleList.isEmpty() ) return null;
//    csVelPickBundle bundlePrev = myCurrentActiveBundleList.get(0);
//    if( myCurrentActiveBundleList.size() > 1 ) {
//      for( int i = 1; i < myCurrentActiveBundleList.size(); i++ ) {
//        csVelPickBundle bundle = myCurrentActiveBundleList.get(i);
//        if( bundle.indexOfTrace1 > traceIndex ) break;
//        bundlePrev = bundle;
//      }
//    }
//    return( bundlePrev.velFunc );
//  }
//  private csEnsemble traceToEns( float traceIndex ) {
//    if( myCurrentActiveBundleList.isEmpty() ) return null;
//    csVelPickBundle bundlePrev = myCurrentActiveBundleList.get(0);
//    if( myCurrentActiveBundleList.size() > 1 ) {
//      for( int i = 1; i < myCurrentActiveBundleList.size(); i++ ) {
//        csVelPickBundle bundle = myCurrentActiveBundleList.get(i);
//        if( bundle.indexOfTrace1 > traceIndex ) break;
//        bundlePrev = bundle;
//      }
//    }
//    return( bundlePrev.ens );
//  }
  @Override
  public synchronized void draw( csSeisView seisview, Graphics2D g ) {
    if( seisview.getTraceBuffer() == null ) return;
    if( myIsBundleLock ) {
      return;
    }
    myIsBundleLock = true;

    RenderingHints rhints_save = g.getRenderingHints();
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

    //-------------------------------------------------------------------
    // Draw velocity functions
    //
    g.setStroke( new BasicStroke(2.0f) );
    // Loop over all currently displayed ensembles:
    for (csVelPickBundle bundle : myCurrentActiveBundleList) {
      for( int ivelField = 0; ivelField < myVelFieldList.size(); ivelField++ ) {
        if( myIsActive && ivelField == myActiveVelFieldIndex ) continue;
        TreeMap<csEnsemble,csVelFunction> map = myVelFieldList.get(ivelField);
        if( !map.containsKey( bundle.ens ) ) continue;
        drawInternal( map.get(bundle.ens).getPicks(), false, g, myDisplayAttrList.get(ivelField), bundle.ens, seisview );
      }
      if( myIsActive && myActiveVelFieldIndex >= 0 ) {
        java.util.List<csVelFunction.Item> list = bundle.velFunc.getPicks();
        drawInternal( list, true, g, myDisplayAttrList.get(myActiveVelFieldIndex), bundle.ens, seisview );
      }
    }
    g.setRenderingHints(rhints_save);
    myIsBundleLock = false;
  }
  private void drawInternal( java.util.List<csVelFunction.Item> picks, boolean plotPoints, Graphics2D g, csHorizonAttr attr, csEnsemble ens, csSeisView seisview ) {
    if( picks.isEmpty() ) return;
    g.setColor( attr.color );

    GeneralPath path = new GeneralPath();
    // Extrapolate from top of gather
    path.moveTo( seisview.xModel2View( velToTrace( picks.get(0).vel,ens) ), seisview.yModel2View( 0 ) );
    for (csVelFunction.Item item : picks) {
      float xView = seisview.xModel2View( velToTrace(item.vel,ens) );
      float yView = seisview.yModel2View( item.time/mySampleInt );
//      System.out.println("Pick " + ens + " "+ item.vel + " " + item.time/mySampleInt + "  " + xView + " " + yView );
      path.lineTo( xView, yView );
      if( plotPoints ) {
        Shape shape = new Ellipse2D.Float( xView-5, yView-5, 11, 11);
        g.draw(shape);
      }
    }
    // Extrapolate to bottom of gather
    path.lineTo( seisview.xModel2View( velToTrace( picks.get(picks.size()-1).vel,ens) ), seisview.yModel2View( myNumSamples-1 ) );
    g.draw( path );
  }
  //-----------------------------------------------------------
  // Mouse and picking routines
  //
  private void pickMouseEvent( MouseEvent e ) {
    Point pointView = e.getPoint();
    csSampleInfo sInfo = mySeisView.getSampleInfo( pointView.x, pointView.y );
    if( sInfo.trace < 0 || sInfo.trace >= mySeisView.getTraceBuffer().numTraces() ) return;
//    System.out.println("pick trace/time " + sInfo.traceDouble + " / " + sInfo.time);
    csVelPickBundle bundle = traceToBundle( (float)sInfo.traceDouble );
    if( bundle == null ) return;
    if( SwingUtilities.isLeftMouseButton(e) ) {
      float vel = bundle.traceToVel( (float)sInfo.traceDouble );
      bundle.velFunc.addPick( (float)sInfo.time*1000.0f, vel );
      mySeisView.repaint();
    }
    else if( SwingUtilities.isMiddleMouseButton(e) ) {
      bundle.velFunc.removePick( (float)sInfo.time*1000.0f );
      mySeisView.repaint();
    }
  }
  @Override
  public void mouseDragged( MouseEvent e ) {
  }
  @Override
  public void mouseClicked(MouseEvent e) {
    if( !myIsActive || myActiveVelFieldIndex < 0 ) return;
    pickMouseEvent( e );
  }
  @Override
  public void mouseReleased(MouseEvent e) {
  }

  @Override
  public void mouseMoved(MouseEvent e) {
  }
  @Override
  public void mousePressed(MouseEvent e) {
  }
  @Override
  public void mouseEntered(MouseEvent e) {
  }
  @Override
  public void mouseExited(MouseEvent e) {
  }

  @Override
  public void keyPressed( KeyEvent e ) {
  }
  @Override
  public void keyReleased( KeyEvent e ) {
  }

}
