/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.seaview;

import cseis.general.csStandard;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;

/**
 * Dialog window showing license information
 * @author 2013 Felipe Punto
 */
@SuppressWarnings("serial")
public class csSeaViewLicenseDialog extends JDialog {
  private static final String DUMMY_ADDRESS = "name@site.com";
  private JButton myButtonYes;
  private JButton myButtonNo;
  private JCheckBox myBoxSendEmail;
  private JTextField myTextEmailAddress;
  private Color myColorEmail;

  private boolean myIsEmailEnabled;
  private boolean myIsAccepted;
  private boolean myIsAcceptedSendEmail;
  private String myEmailAddress;
  
  public csSeaViewLicenseDialog() {
    setTitle( "License information" );
    setModal( true );
    setAlwaysOnTop(true);

    myColorEmail = new Color(255,200,100);
    
    String date    = "15 July 2013";
    String message = "Running the " + date + " version (Release " + SeaView.VERSION + ") of OpenSeaSeis.";
    String address = "john@dix.mines.edu";
    String textEmail =
            "To give us at CWP an idea of who uses OpenSeaSeis the following message\n" +
            "     '" + message + "'\n" +
            "will be emailed to: " + address + "\n" +
            "You will then be put on the OpenSeaSeis e-mail list" +
            "and will be informed of future updates of OpenSeaSeis. " +
            "However, if you would rather not have this message sent " +
            "you may unselect the check box below.";
    myIsAccepted = true;
    myIsAcceptedSendEmail = true;
    myBoxSendEmail = new JCheckBox("Send automatic mail message to CWP");
    myBoxSendEmail.setSelected(true);
    JTextArea textAreaEmail = new JTextArea( textEmail );
    textAreaEmail.setEditable(false);
    textAreaEmail.setLineWrap(true);
    textAreaEmail.setWrapStyleWord(true);
    Font font = textAreaEmail.getFont();
    textAreaEmail.setFont( new Font(font.getName(),Font.ITALIC,font.getSize()) );
    myTextEmailAddress = new JTextField( DUMMY_ADDRESS );
    myTextEmailAddress.setToolTipText("Specify your email address here");
    myTextEmailAddress.setBackground( myColorEmail );
    myTextEmailAddress.setBorder( BorderFactory.createLineBorder(Color.black) );
    myTextEmailAddress.setFocusable(true);
    
    JLabel labelAddress = new JLabel("Specify your email address here:  ");
    JLabel labelAccept = new JLabel(
            "<html><i><b>By answering you agree to abide by the terms and conditions<br>"
                + "of the above LEGAL STATEMENT:</b?</i></html>", JLabel.CENTER);
    
    myButtonYes = new JButton("Yes");
    myButtonNo  = new JButton("No");
    
    JTextArea area = new JTextArea();
    area.setFont( new Font(font.getName(),Font.PLAIN,font.getSize()-1) );
    area.setEditable(false);
    area.setLineWrap(false);
    String text = cseis.resources.csResources.getText( "OpenSeaSeis_LEGAL_STATEMENT" );
    if( text == null ) {
      area.append("*** COULD NOT FIND LICENSE INFORMATION FILE ***");
    }
    else {
      area.append(text);
    }
    JScrollPane pane = new JScrollPane( area,
            JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED );

    JPanel panelLicense = new JPanel( new BorderLayout() );
    panelLicense.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("LEGAL STATEMENT"),
        csStandard.INNER_EMPTY_BORDER ) );

    JPanel panelEmail = new JPanel( new BorderLayout() );
    panelEmail.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createTitledBorder("EMAIL"),
        csStandard.INNER_EMPTY_BORDER ) );
    JPanel panelSendEmail = new JPanel( new BorderLayout() );

    panelLicense.add(pane,BorderLayout.CENTER);
    
    panelSendEmail.add(myBoxSendEmail, BorderLayout.NORTH);
    panelSendEmail.add(labelAddress, BorderLayout.WEST);
    panelSendEmail.add(myTextEmailAddress, BorderLayout.CENTER);
    panelEmail.add(panelSendEmail, BorderLayout.NORTH);
    panelEmail.add(textAreaEmail, BorderLayout.CENTER);
    
    JPanel panelButtons = new JPanel( new GridBagLayout() );
    int xp = 0;
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( labelAccept, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.HORIZONTAL, new Insets( 11, 10, 0, 20 ), 0, 0 ) );
    panelButtons.add( myButtonYes, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTH,
        GridBagConstraints.NONE, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( myButtonNo, new GridBagConstraints(
        xp++, 0, 1, 1, 0.0, 0.0, GridBagConstraints.SOUTHWEST,
        GridBagConstraints.NONE, new Insets( 11, 0, 0, 0 ), 0, 0 ) );
    panelButtons.add( Box.createHorizontalGlue(), new GridBagConstraints(
        xp++, 0, 1, 1, 0.5, 1.0, GridBagConstraints.CENTER,
        GridBagConstraints.HORIZONTAL, new Insets( 0, 0, 0, 0 ), 0, 0 ) );

    JPanel panelTop = new JPanel(new BorderLayout());
    panelTop.add( panelLicense, BorderLayout.CENTER );
    // Enable the following line to include email in dialog:
    myIsEmailEnabled = false;
    //    panelTop.add( panelEmail, BorderLayout.SOUTH );

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( csStandard.DIALOG_BORDER );

    panelAll.add(panelTop,BorderLayout.CENTER);
    panelAll.add(panelButtons,BorderLayout.SOUTH);
    getContentPane().add(panelAll);
    
    Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
    double width = screenSize.getWidth();
    double height = screenSize.getHeight();
    setSize( new Dimension(600,600) );
    setLocation( (int)((width-600)/2), (int)((height-600)/2) );
    
    myButtonYes.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        if( apply() ) {
          myIsAccepted = true;
          cancel();
        }
      }
    });
    myButtonNo.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        myIsAccepted = false;
        cancel();
      }
    });
    myBoxSendEmail.addActionListener( new ActionListener() {
      @Override
      public void actionPerformed( ActionEvent e ) {
        boolean isSelected = myBoxSendEmail.isSelected();
        myTextEmailAddress.setEnabled( isSelected );
        if( isSelected ) {
          myTextEmailAddress.setBackground( myColorEmail );
          myTextEmailAddress.setBorder( BorderFactory.createLineBorder(Color.black) );
        }
        else {
          myTextEmailAddress.setBackground( myButtonYes.getBackground() );
          myTextEmailAddress.setBorder( null );
        }
      }
    });
    addWindowListener(new WindowAdapter(){ 
      @Override
      public void windowOpened( WindowEvent e){ 
        myTextEmailAddress.requestFocusInWindow();
      } 
    });     
  }
  private boolean apply() {
    if( !myIsEmailEnabled ) return true;
    myIsAcceptedSendEmail = myBoxSendEmail.isSelected();
    if( myIsAcceptedSendEmail ) {
      myEmailAddress = myTextEmailAddress.getText();
      if( myEmailAddress.length() < 3 || !myEmailAddress.contains("@") || myEmailAddress.compareTo(DUMMY_ADDRESS) == 0 ) {
        setAlwaysOnTop(false);
        String message = "";
        if( myEmailAddress.length() == 0 ) {
          message += "No email address specified.\n";
        }
        else {
          message += "Specified email address does not appear to be valid:\n" +
                  "   '" + myEmailAddress + "'\n";
        }
        message += "Please provide your email address in the text field.";
        JOptionPane.showMessageDialog( this, message,
                "Error", JOptionPane.ERROR_MESSAGE );
        setAlwaysOnTop(true);
        return false;
      }
    }
    return true;
  }
  private void cancel() {
    dispose();
  }
  public boolean acceptedLicense() {
    return myIsAccepted;
  }
  public boolean acceptedSendEmail() {
    return myIsAcceptedSendEmail;
  }
  public String getEmailAddress() {
    return myEmailAddress;
  }
  public static void main( String[] args ) {
    csSeaViewLicenseDialog dialog = new csSeaViewLicenseDialog();
    dialog.setVisible(true);
  }
}

