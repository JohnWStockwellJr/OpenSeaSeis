/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

package cseis.xcseis;

import cseis.jni.csNativeModuleHelp;
import cseis.seaview.SeaView;
import cseis.swing.csFileMenuEvent;
import cseis.swing.csFileMenuListener;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.*;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 *
 * @author bjorn
 */
public class XCSeis extends JFrame implements csFileMenuListener, csIFlowViewListener {
  public static final String VERSION = "0.20";
  public static final int DEFAULT_DIVIDER_SIZE = 10;
  private static final String TEXT_NEW_FLOW = "NewFlow";
  private JButton myButtonFormat;
  private JButton myButtonSubmit;
  private JToggleButton myButtonFlowView;
  private csPreferences myPreferences;
  private Font mySystemFont;
  
  // Flows
  private csFlowPanel myActivePanel;
//  private String myPathActiveFlow;
//  private ArrayList<csFlowBundle> myFlows;
  private csProjectDef myProjDef;
//  private ArrayList<csFlowPanel> myFlowPanelList;
  private JTabbedPane myTabPaneFlows;
  private int myTabIndexModuleHelpPane;

  // Menu & status bar
  private csMenuBar myMenuBar;

  private JToolBar myToolBar;
  private csStatusBar myStatusBar;
  private JFileChooser myFileChooserFlows;

  // Special text areas
//  private JTextPane myTextPaneFlow;
  private JTextArea myTextAreaLog;
  private JTextArea myTextAreaParam;
  private JTextArea myTextOutput;

  // Panels
  private final csSubPanel myPanelFlow;
  private final csSubPanel myPanelLog;
  private final JPanel myPanelParam;
  private final JPanel myPanelHelp;
  private final JPanel myPanelModuleList;
  private final JPanel myPanelHeaderList;
  private final csSubPanel myPanelOutput;

  // Scroll panes
  private JScrollPane myPaneLog;
  private JScrollPane myPaneOutput;
  private csModuleParamPanel myPaneParam;
  private csModuleHelpPanel myPaneHelp;
  private csModuleListingPanel myPaneModuleList;
  private csStandardHeaderPanel myPaneHeaderList;
  private JScrollPane myPanePreview;

  // Combined panels:
  private JSplitPane mySplitPaneParamHelp;
  private JSplitPane mySplitPaneFlowParamHelp;
  private JSplitPane mySplitPaneLogOutput;
  private JSplitPane mySplitPaneAll;
  private JTabbedPane myTabPaneHelp;
  private JTabbedPane myTabPaneModuleParam;

  public XCSeis() {
    super("XCSeis, version " + XCSeis.VERSION);
    
    mySystemFont = new JLabel().getFont();
    myPreferences = new csPreferences( this );
    myPreferences.readPreferences();
//    myFlowPanelList = new ArrayList<csFlowPanel> (1);

    myStatusBar = new csStatusBar();

    myMenuBar = new csMenuBar(this);
    java.util.List list = myPreferences.getRecentFileList();
    if( list != null ) {
      for( int i = list.size()-1; i >= 0; i-- ) {  // Add files backwards. Last added file goes to top of list
        myMenuBar.addRecentFile( (String)list.get(i) );
      }
    }
    myProjDef = new csProjectDef( myPreferences.getProjDirectory() );
    myFileChooserFlows = new JFileChooser( myProjDef.pathFlows() );
  
    myTextOutput = new JTextArea();
    myTextOutput.setEditable(false);
    myTextOutput.setWrapStyleWord(true);
    myTextOutput.setWrapStyleWord(true);
//    try {
//      myTextOutput.read( new FileInputStream(myFileStdout), "");
//    }
//    catch (IOException ex) {
//      myTextOutput.setText("...problem occurred when capturing standard output stream.");
//    }

    myButtonFlowView = new JToggleButton("Flow View     ", false);
    myButtonSubmit = new JButton("Submit Flow   ");
    
    myTextAreaLog   = new JTextArea();
    myTextAreaParam = new JTextArea();

    myTextAreaLog.setEditable( false );
    myTextAreaLog.setLineWrap( false );
    myTextAreaLog.setFont( new Font("Monospaced",Font.PLAIN, 11 ));
    myTextAreaParam.setEditable( true );
    myTextAreaParam.setLineWrap( false );

    myPanelFlow  = new csSubPanel( "<html> &nbsp; <b>Flow editor</b> </html>" );
    myPanelLog   = new csSubPanel( "<html> &nbsp; <b>Log viewer</b></html>" );
    myPanelOutput = new csSubPanel( "<html> &nbsp; <b>Output</b></html>" );
//    myPanePreview = new csSubPanel( "<html> &nbsp; <b>Preview</b></html>" );
    myPanelParam = new JPanel( new BorderLayout() );
    myPanelHelp  = new JPanel( new BorderLayout() );
    myPanelModuleList = new JPanel( new BorderLayout() );
    myPanelHeaderList = new JPanel( new BorderLayout() );
//    myPanelPreview = new JPanel( new BorderLayout() );

    csFlowPanel panel = new csFlowPanel(this,TEXT_NEW_FLOW);
    panel.updateModel(new StringBuffer(""));
    myTabPaneFlows  = new JTabbedPane();
    myTabPaneFlows.add( TEXT_NEW_FLOW, panel );
    myTabPaneFlows.setTabComponentAt(myTabPaneFlows.getTabCount()-1,new csButtonTab(this,myTabPaneFlows));
    myTabPaneFlows.addChangeListener( new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
        int selectedIndex = myTabPaneFlows.getSelectedIndex();
        updateActiveFlow(  );
      }
    });
    
    myTabPaneModuleParam = new JTabbedPane( JTabbedPane.TOP, JTabbedPane.SCROLL_TAB_LAYOUT );
    myTabPaneModuleParam.addTab( "Parameter editor", myPanelParam );
    myTabPaneModuleParam.addTab( "Module list", myPanelModuleList );
    myTabPaneModuleParam.setSelectedIndex( myTabPaneModuleParam.getTabCount()-1 );

    myTabPaneHelp = new JTabbedPane( JTabbedPane.TOP, JTabbedPane.SCROLL_TAB_LAYOUT );
    myTabIndexModuleHelpPane = 0;
    myTabPaneHelp.addTab( "Help", myPanelHelp );
    myTabPaneHelp.addTab( "Header list", myPanelHeaderList );

    // Scroll panes
    myPaneLog        = new JScrollPane( myTextAreaLog, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );
    myPaneOutput     = new JScrollPane( myTextOutput, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED );
    myPaneParam      = new csModuleParamPanel(); //( myTextAreaParam, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );
    myPaneHelp       = new csModuleHelpPanel(); //( myLabelHelp, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );
    myPaneModuleList = new csModuleListingPanel(); //( myLabelModuleList, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );
    myPaneHeaderList = new csStandardHeaderPanel(); //( myLabelModuleList, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS );
    myPanePreview    = new JScrollPane(  );

    // Combined panels:
    myPanelFlow.add( myTabPaneFlows,  BorderLayout.CENTER );
    myPanelLog.add( myPaneLog,  BorderLayout.CENTER );
    myPanelOutput.add( myPaneOutput,  BorderLayout.CENTER );
//    myPanelParam.add( myTitleParam, BorderLayout.NORTH );
    myPanelParam.add( myPaneParam,  BorderLayout.CENTER );
    myPanelHelp.add( myPaneHelp,  BorderLayout.CENTER );
    myPanelModuleList.add( myPaneModuleList,  BorderLayout.CENTER );
    myPanelHeaderList.add( myPaneHeaderList,  BorderLayout.CENTER );
//    myPanelPreview.add( myPanePreview,  BorderLayout.CENTER );

    mySplitPaneParamHelp     = new JSplitPane( JSplitPane.VERTICAL_SPLIT, myTabPaneModuleParam, myTabPaneHelp );
    mySplitPaneFlowParamHelp = new JSplitPane( JSplitPane.HORIZONTAL_SPLIT, mySplitPaneParamHelp, myPanelFlow );
    mySplitPaneLogOutput = new JSplitPane( JSplitPane.VERTICAL_SPLIT, myPanelLog, myPanelOutput );
    mySplitPaneAll      = new JSplitPane( JSplitPane.HORIZONTAL_SPLIT, mySplitPaneFlowParamHelp, mySplitPaneLogOutput );

    mySplitPaneParamHelp.setDividerSize( DEFAULT_DIVIDER_SIZE );
    mySplitPaneFlowParamHelp.setDividerSize( DEFAULT_DIVIDER_SIZE );
    mySplitPaneAll.setDividerSize( DEFAULT_DIVIDER_SIZE );
    mySplitPaneLogOutput.setDividerSize( DEFAULT_DIVIDER_SIZE );

    mySplitPaneParamHelp.setOneTouchExpandable(true);
    mySplitPaneFlowParamHelp.setOneTouchExpandable(true);
    mySplitPaneLogOutput.setOneTouchExpandable(true);
    mySplitPaneAll.setOneTouchExpandable(true);

    myButtonFormat = new JButton("Format");
    myToolBar = new JToolBar();
//    myToolBar.add(myButtonFormat);
//    myToolBar.addSeparator();
    myToolBar.add(myButtonFlowView);
    myToolBar.addSeparator();
    myToolBar.add(myButtonSubmit);

    setJMenuBar( myMenuBar );
    getContentPane().add( myToolBar, BorderLayout.NORTH );
    getContentPane().add( mySplitPaneAll, BorderLayout.CENTER );
    getContentPane().add( myStatusBar , BorderLayout.SOUTH );

//    ImageIcon appIcon = cseis.resources.csResources.getIcon("XCSeis.png");
//    if( appIcon != null ) setIconImage( appIcon.getImage() );

    pack();
    Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    Dimension appSize = new Dimension( screenSize.width * 9 / 10, screenSize.height * 8 / 10 );
    setSize( appSize );
    setLocation( ( screenSize.width - appSize.width ) / 2, ( screenSize.height - appSize.height ) / 2 );
    setSize( appSize );
    setLocation( ( screenSize.width - appSize.width ) / 2, 0 );
    mySplitPaneParamHelp.setDividerLocation( (int)(appSize.height / 3) );
    mySplitPaneFlowParamHelp.setDividerLocation( (int)(appSize.width / 4) );
    mySplitPaneLogOutput.setDividerLocation( (int)(2*appSize.height / 3) );
    mySplitPaneAll.setDividerLocation( (int)(3*appSize.width / 4) );

    myPaneModuleList.addModuleListingListener( new csIModuleListingListener() {
      public void selectModule(String moduleName) {
        showModuleHelp( moduleName );
      }
      public void addModule( String moduleName ) {
        csNativeModuleHelp moduleHelp = new csNativeModuleHelp();
        String text = moduleHelp.moduleExample(moduleName) + "\n"; // Add empty line at end
        if( myActivePanel == null ) {
          updateActiveFlow();
        }
        if( myActivePanel != null ) {
          myActivePanel.addModule( moduleName, text );
        }
      }
    });
      
    myButtonFormat.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        // myFlowPanelList.get(myTabPaneFlows.getSeulelectedIndex()).formatText(true);
      }
    });
    myButtonSubmit.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        if( myActivePanel != null ) {
          submitFlow();
        }
      }
    });
    myButtonFlowView.addActionListener( new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        boolean selected = myButtonFlowView.isSelected();
        for( int i = 0; i < myTabPaneFlows.getTabCount(); i++ ) {
          csFlowPanel panel = (csFlowPanel)myTabPaneFlows.getComponentAt(i);
          panel.setFlowView(selected);
        }
      }
    });

    super.addWindowListener(new WindowAdapter() {
      public void windowClosing( WindowEvent e ) {
        exitApplication();
      }
    });
  }
  @Override
  public void setVisible( boolean doSet ) {
    super.setVisible(doSet);
    if( !myProjDef.checkDirs() ) {
      projectSetup();
    }
  }
  public boolean refreshFile() {
    if( myActivePanel != null ) {
      String filename = myActivePanel.flowPathName();
      closeFlow();
      return openFlow( filename );
    }
    return false;
  }
  public void openSettingsDialog() {
    csSettingsDialog dialog = new csSettingsDialog(this);
    dialog.setVisible(true);
  }
  public Font getSystemFont() {
    return mySystemFont;
  }
  public void setSystemFont( Font font ) {
    mySystemFont = font;
    if( myTabPaneFlows != null ) {
      for( int itab = 0; itab < myTabPaneFlows.getTabCount(); itab++ ) {
        csFlowPanel p = (csFlowPanel)myTabPaneFlows.getComponentAt(itab);
        p.updateFont(font);
      }
    }
  }
  public String pathActiveFlow() {
    return myActivePanel.flowPathName();
  }
  public boolean closeFlow() {
    return closeFlow( myTabPaneFlows.getSelectedIndex() );
  }
  public boolean closeFlow( int selectedIndex ) {
    if( selectedIndex >= 0 && selectedIndex < myTabPaneFlows.getTabCount() ) {
     csFlowPanel p = (csFlowPanel)myTabPaneFlows.getComponentAt(selectedIndex);
     //csFlowPanel p = myFlowPanelList.get(selectedIndex);
      if( p.isUpdated() ) {
        int option = JOptionPane.showConfirmDialog(this, "Save file before closing?",
                "Confirm dialog",JOptionPane.YES_NO_CANCEL_OPTION);
        if( option == JOptionPane.YES_OPTION ) {
          myTabPaneFlows.setSelectedIndex(selectedIndex);
          if( p.getText().compareTo('*'+TEXT_NEW_FLOW) == 0 ) {
            saveFlowAs();
          }
          else {
            saveFlow();
          }
        }
        else if( option == JOptionPane.NO_OPTION ) {
//          myTabPaneFlows.setTitleAt(selectedIndex,stripFilename(myPathActiveFlow));
        }
        else if( option == JOptionPane.CANCEL_OPTION ) {
          return false;
        }
      }
      myTabPaneFlows.remove(selectedIndex);
//      myFlowPanelList.remove(selectedIndex);
      if( selectedIndex >= myTabPaneFlows.getTabCount() ) selectedIndex -= 1;
      if( selectedIndex >= 0 ) {
        myTabPaneFlows.setSelectedIndex( selectedIndex );
        updateActiveFlow();
      }
      return true;
    }
    return false;
  }
  public boolean openFlow() {
    JFileChooser fc = myFileChooserFlows;
    int option = fc.showOpenDialog( XCSeis.this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      String filename = fc.getSelectedFile().getAbsolutePath();
      fc.setCurrentDirectory( fc.getCurrentDirectory() );
      return openFlow(filename);
    }
    return false;
  }
  public void newFlow() {
    csFlowPanel panel = new csFlowPanel(this,TEXT_NEW_FLOW);
    panel.updateModel(new StringBuffer(""));
//    myFlowPanelList.add(panel);
    myTabPaneFlows.add( TEXT_NEW_FLOW, panel );
    myTabPaneFlows.setTabComponentAt(myTabPaneFlows.getTabCount()-1,new csButtonTab(this,myTabPaneFlows));
    updateActiveFlow();
  }
  public void updateStatus(String text) {
    myStatusBar.setMessage(text);
  }
  public void updateFlowTab( csFlowPanel panelIn ) {
    for( int itab = 0; itab < myTabPaneFlows.getTabCount(); itab++ ) {
      csFlowPanel p = (csFlowPanel)myTabPaneFlows.getComponentAt(itab);
      if( p.equals(panelIn) ) {
        myTabPaneFlows.setTitleAt( itab, '*'+myTabPaneFlows.getTitleAt(itab) );
        myTabPaneFlows.repaint();
        updateActiveFlow();
        break;
      }
    }
  }
  public boolean openFlow( String filename ) {
    // First, check if file is already open. If it is, switch to that tab and refresh
    for( int itab = 0; itab < myTabPaneFlows.getTabCount(); itab++ ) {
      csFlowPanel panel = (csFlowPanel)myTabPaneFlows.getComponentAt(itab);
      if( filename.compareTo(panel.flowPathName()) == 0 ) {
        myTabPaneFlows.setSelectedIndex(itab);
        setTitle(filename);
        return true;
      }
    }
    try {
      if( filename == null ) return false;
      FileReader fileReader = new FileReader(filename);
      if( fileReader == null ) return false;
      BufferedReader reader = new BufferedReader( fileReader );
      String line;
      StringBuffer strBuffer = new StringBuffer();
      while( ( line = reader.readLine() ) != null ) {
        strBuffer.append(line + "\n");
      }
      csFlowPanel panel = new csFlowPanel(this,filename);
      myActivePanel = panel;
      panel.updateModel(strBuffer);
      String title = stripFilename(filename);
      if( myTabPaneFlows.getSelectedIndex() >= 0 && myTabPaneFlows.getTitleAt(myTabPaneFlows.getSelectedIndex()) == this.TEXT_NEW_FLOW ) {
//        myFlowPanelList.set(myTabPaneFlows.getSelectedIndex(),panel);
        myTabPaneFlows.setComponentAt(myTabPaneFlows.getSelectedIndex(), panel);
        myTabPaneFlows.setTitleAt(myTabPaneFlows.getSelectedIndex(),title);
        myTabPaneFlows.setTabComponentAt(myTabPaneFlows.getSelectedIndex(),new csButtonTab(this,myTabPaneFlows));
      }
      else {
//        myFlowPanelList.add(panel);
        myTabPaneFlows.add(title,panel);
        myTabPaneFlows.setSelectedIndex(myTabPaneFlows.getTabCount()-1);
        myTabPaneFlows.setTabComponentAt(myTabPaneFlows.getTabCount()-1,new csButtonTab(this,myTabPaneFlows));
      }
      if( myButtonFlowView.isSelected() ) panel.setFlowView(true);
      myMenuBar.addRecentFile( filename );
    }
    catch( IOException exc ) {
      exc.printStackTrace();
      return false;
    }
    updateActiveFlow();
    return true;
  }
  private void updateActiveFlow() {
    int selectedIndex = myTabPaneFlows.getSelectedIndex();
    int enableFlag = csMenuBar.NONE;
    if( selectedIndex >= 0 ) {
      csFlowPanel panel = (csFlowPanel)myTabPaneFlows.getSelectedComponent();
      myActivePanel = panel;
      String text = myTabPaneFlows.getTitleAt(selectedIndex);
      if( text.compareTo(TEXT_NEW_FLOW) == 0 ) {
        enableFlag = csMenuBar.NEW_FILE;
      }
      else if( text.compareTo('*'+TEXT_NEW_FLOW) == 0 ) {
        enableFlag = csMenuBar.NEW_CHANGED_FILE;
      }
      else if( text.length() > 0 && text.charAt(0) == '*' ) {
        enableFlag = csMenuBar.EXISTING_CHANGED_FILE;
      }
      else {
        enableFlag = csMenuBar.EXISTING_FILE;
      }
    }
    else {
      myActivePanel = null;
    }
    myMenuBar.setMenusEnabled( enableFlag );
    if( myActivePanel != null ) {
      setTitle( myActivePanel.flowPathName() );
    }
    else {
      setTitle("");
    }
  }
  private String stripFilename( String pathName ) {
    int counter = pathName.length() - 1;
    while( counter >= 0 && pathName.charAt(counter) != '/') {
      counter -= 1;
    }
    return pathName.substring(counter+1);
  }
  public void saveFlow() {
    if( myActivePanel == null ) return;
    BufferedWriter writer = null;
    try {
      writer = new BufferedWriter( new FileWriter( myActivePanel.flowPathName() ) );
      csFlowPanel panel = (csFlowPanel)myTabPaneFlows.getSelectedComponent();
      if( panel == null ) return;
      writer.write( panel.getText() );
      writer.close();
      updateStatus("...saved flow " + myActivePanel.flowPathName());
      myTabPaneFlows.setTitleAt(myTabPaneFlows.getSelectedIndex(),stripFilename(myActivePanel.flowPathName()));
      ((csFlowPanel)myTabPaneFlows.getSelectedComponent()).resetUpdated();
      updateActiveFlow();
    }
    catch( IOException e ) {
      e.printStackTrace();
    }
  }
  public void saveFlowAs() {
    JFileChooser fc = myFileChooserFlows;
    int option = fc.showSaveDialog( XCSeis.this );
    if( option == JFileChooser.APPROVE_OPTION ) {
      String filename = fc.getSelectedFile().getAbsolutePath();
      fc.setCurrentDirectory( fc.getCurrentDirectory() );
      setTitle( filename );
      try {
        BufferedWriter writer = new BufferedWriter( new FileWriter( filename ) );
        csFlowPanel panel = (csFlowPanel)myTabPaneFlows.getSelectedComponent();
        myActivePanel = panel;
        if( panel == null ) return;
        writer.write( panel.getText() );
        writer.close();
        myTabPaneFlows.setTitleAt(myTabPaneFlows.getSelectedIndex(),stripFilename(filename));
        ((csFlowPanel)myTabPaneFlows.getSelectedComponent()).resetUpdated();
        updateStatus("...saved flow as " + filename);
        myMenuBar.addRecentFile( filename );
        updateActiveFlow();
      }
      catch( IOException exc ) {
        exc.printStackTrace();
        return;
      }
    }
  }
  public void openInSeaView( final String filename ) {  
    String path = System.getProperty( "java.library.path" );
    String libName = System.mapLibraryName( "csJNIlib" );
    try {
      System.load( path + "/" + libName );
    }
    catch( java.lang.UnsatisfiedLinkError e ) {
      JOptionPane.showMessageDialog( null,
          e.toString() + "\n" +
          "java.library.path = " + System.getProperty( "java.library.path" ) + "\n" +
          " - Seaview will not run.", "Error",
          JOptionPane.ERROR_MESSAGE );
      return;
    }

    final SeaView seaview = new SeaView( true );
    seaview.setVisible(true);

    SwingUtilities.invokeLater( new Runnable() {
      public void run() {
        seaview.openFile(filename);
      }
    });
  }
  public void projectSetup() {
    final csProjectDefDialog dialog = new csProjectDefDialog( this, myProjDef );
    dialog.setVisible(true);
  }
  public void submitFlow() {
    myPanelLog.setInfo("..running..", csSubPanel.STATUS_PENDING );
    updateStatus("...submitting flow " + myActivePanel);
    Thread thread = new Thread() {
      @Override
      public void run() {
        String flowPath = myActivePanel.flowPathName();
        String flowName = flowPath;
        String flowRootPath = null;
        int index = flowPath.lastIndexOf("/");
        if( index >= 0 ) {
          flowName = flowPath.substring(index);
          flowRootPath = flowPath.substring(0,index);
        }
        if( flowRootPath == null ) {
          JOptionPane.showMessageDialog(rootPane, "No flow file: Make sure to save your flow to disk before job submission", "Error", JOptionPane.ERROR_MESSAGE);
          return;
        }
        String logName = myProjDef.pathLogs() + flowName.substring(0,flowName.length()-5) + ".log";
        String commandText = "seaseis -f " + flowPath + " -o " + logName;

        writeLineOut( "Submit command:" );
        writeLineOut( commandText );
        try {
          final Process process = Runtime.getRuntime().exec(commandText, null, new File(flowRootPath) );
          InputStream stderr = process.getErrorStream();
          InputStream stdout = process.getInputStream();

          // clean up if any output in stdout
          BufferedReader brCleanUp = new BufferedReader (new InputStreamReader (stdout));
          String line;
          while( (line = brCleanUp.readLine ()) != null ) {
            writeLineOut(line);
          }
          brCleanUp.close();
          brCleanUp = new BufferedReader (new InputStreamReader (stderr));
          while( (line = brCleanUp.readLine()) != null ) {
            writeLineOut(line);
          }
          brCleanUp.close();
          try {
            process.waitFor();
            int exitValue = process.exitValue();
            if( exitValue != 0 ) {
              myPanelLog.setInfo("..error..", csSubPanel.STATUS_ERROR );
            }
            else {
              myPanelLog.setInfo("..successful..", csSubPanel.STATUS_OK );
            }
            writeLineOut( "EXIT value: " + exitValue );
          }
          catch (InterruptedException ex) {
          }
        }
        catch( Exception err ) {
          JOptionPane.showMessageDialog(rootPane, err.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
        }
        updateLogFile(logName);
        setSubmissionEnabled( true );
      }
    };
    setSubmissionEnabled( false );
    thread.start();
  }
  private void setSubmissionEnabled( boolean doEnable ) {
    myButtonSubmit.setEnabled(doEnable);
    myMenuBar.setSubmissionEnabled(doEnable);
  }
  private void writeLineOut( String line ) {
    myTextOutput.append( line + "\n" );
  }
  public void openRecentLog() {
    if( myActivePanel == null ) return;
    String flowName = myActivePanel.flowPathName();
    int index = flowName.lastIndexOf("/");
    if( index >= 0 ) {
      flowName = flowName.substring(index);
    }
    String logName = myProjDef.pathLogs() + flowName.substring(0,flowName.length()-5) + ".log";
    updateLogFile(logName);
  }
  public void updateLogFile( String logName ) {
    try {
      BufferedReader reader = new BufferedReader(new FileReader(logName));
      String line = null;
      myTextAreaLog.setText("");
      try {
        while ((line = reader.readLine()) != null) {
          myTextAreaLog.append(line + "\n");
        }
      }
      catch( IOException e2 ) {
      }
    }
    catch( FileNotFoundException e ) {
      JOptionPane.showMessageDialog(this, "Cannot open log file: " + e.getMessage(),
              "Error", JOptionPane.ERROR_MESSAGE);
    }
  }
  //------------------------------------------------
  public void exitApplication() {
    myPreferences.setProjDirectoryPath( myProjDef.pathProj() );
    myPreferences.setRecentFileList(myMenuBar.getRecentFileList());
    myPreferences.writePreferences();
    System.exit(0);
  }

  // Recent file menu listener
  public void fileSelected( csFileMenuEvent e ) {
    File file = e.file();
    boolean success = openFlow( file.getAbsolutePath() );
    if( !success ) {
      myMenuBar.removeRecentFile( file.getAbsolutePath() );
    }
  }
  //---------------------------------------------------------------------------
  public void selectModuleInFlow( int moduleIndex, String moduleName ) {
    showModuleHelp( moduleName );
  }
  public void deleteModuleInFlow( int moduleIndex, String moduleName ) {
    // Not implemented yet
  }
  public void showModuleHelp( String moduleName ) {
    myPaneHelp.updateModule(moduleName);
    csNativeModuleHelp moduleHelp = new csNativeModuleHelp();
    String text = moduleHelp.moduleHtmlHelp(moduleName);
    myPaneHelp.updateModuleHelp(text);
    myPaneParam.updateModuleExample(moduleName);
    myTabPaneHelp.setSelectedIndex( myTabIndexModuleHelpPane );
  }
  public boolean setProjectDir( String projDir ) {
    try {
      myProjDef.set(projDir);
      myProjDef.createDirs();
      myFileChooserFlows.setCurrentDirectory( new File(myProjDef.pathFlows()) );
    }
    catch (Exception ex) {
      JOptionPane.showMessageDialog(this, ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
      return false;
    }
    return true;
  }

  //*********************************************************************************
  //*********************************************************************************
  //*********************************************************************************
  //
  public static void main( String[] args ) {
    String path = System.getProperty( "java.library.path" );
    String libName1 = System.mapLibraryName( "XCSeislib" );
    try {
      System.load( path + "/" + libName1 );
    }
    catch( java.lang.UnsatisfiedLinkError e ) {
      JOptionPane.showMessageDialog( null,
                                     e.toString() + "\n" +
                                     "java.library.path = " + System.getProperty( "java.library.path" ) + "\n" +
                                     " - XCSeis will not show any module documentation.", "Error",
                                     JOptionPane.ERROR_MESSAGE );
    }

    try {
      UIManager.setLookAndFeel( "com.sun.java.swing.plaf.windows.WindowsLookAndFeel" );
    }
    catch( Exception e ) {
      // Nothing
    }
    
    boolean createdNewDirectory = false;
    try {
      createdNewDirectory = cseis.general.csAbstractPreference.setCseisDirectory();
    }
    catch( Exception e ) {
      JOptionPane.showMessageDialog( null,
          e.toString() + " - XCSeis will not run.", "Error",
          JOptionPane.ERROR_MESSAGE );
      System.exit( -1 );
    }
    XCSeis xcseis = new XCSeis();
    xcseis.setVisible(true);
    if( createdNewDirectory ) {
      JOptionPane.showMessageDialog( xcseis,
          "The directory\n" + "'" + cseis.general.csAbstractPreference.getCseisDirectory() + "'\n" +
          "has been created on your system.\n\n" +
          "This directory is used to store configuration files.",
          "Info", JOptionPane.INFORMATION_MESSAGE );
    }
  }
}


