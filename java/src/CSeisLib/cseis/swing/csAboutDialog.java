/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */


package cseis.swing;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridLayout;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JTextArea;

import cseis.general.csStandard;

/**
 * Simple "About" dialog window
 * @author Bjorn Olofsson
 */
@SuppressWarnings("serial")
public class csAboutDialog extends JDialog {
  private String myVersion;
  private String myAuthor;
  private String myProductName;
  private String myDate;
  private String myCompany;
  private String myContact;
  private String myComments;
  private boolean myIsInitialized;
  private ImageIcon myLogo;
  private JFrame myParentFrame;
  
  public csAboutDialog( JFrame parentFrame, String productName, String versionString, String author, String company, String dateString, ImageIcon logo ) {
    this( parentFrame, "About" );
    myParentFrame = parentFrame;
    myProductName = productName;
    myLogo    = logo;
    myVersion = versionString;
    myAuthor  = author;
    myDate    = dateString;
    myCompany = company;
  }
  public csAboutDialog( JFrame parentFrame, String productName ) {
    super( parentFrame, "About" );
    myParentFrame = parentFrame;
    myProductName = productName;
    myVersion     = "";
    myAuthor      = "";
    myDate        = "";
    myCompany     = "";
    myContact     = "";
    myComments    = "";
    myIsInitialized = false;
    myLogo        = null;
  }
  public void setVisible( boolean doSet ) {
    if( !myIsInitialized ) {
      initialize();
    }
    super.setVisible( doSet );
  }
  private void initialize() {
    myIsInitialized = true;
    
    JLabel label1 = new JLabel( "   " + myProductName );
    label1.setFont( new Font( "SansSerif", Font.BOLD, 12 ) );
    JLabel label2 = new JLabel( "Version:   " + myVersion, JLabel.LEFT );
    JLabel label3 = new JLabel( "Author:    " + myAuthor, JLabel.LEFT );

    JPanel panelIcon = new JPanel( new BorderLayout() );
    panelIcon.setBorder( csStandard.INNER_EMPTY_BORDER );
    if( myLogo != null ) {
      JLabel label = new JLabel( myLogo, JLabel.LEFT );
      panelIcon.add( label, BorderLayout.WEST );
    }
    panelIcon.add( label1, BorderLayout.CENTER );
    JPanel panel = new JPanel(new GridLayout(4,1));
    panel.setBorder( csStandard.INNER_EMPTY_BORDER );
    panel.add(label2);
    if( myAuthor.length() != 0 ) panel.add(label3);
    if( myCompany.length() != 0 ) panel.add(new JLabel( "Company:   " + myCompany, JLabel.LEFT ));
    if( myContact.length() != 0 ) panel.add(new JLabel( "Contact:   " + myContact, JLabel.LEFT ));
    if( myDate.length() != 0 )    panel.add(new JLabel( "Date:      " + myDate, JLabel.LEFT ));
    JPanel panel_comments = new JPanel();
    if( myComments.length() != 0 ) {
      JTextArea area = new JTextArea( myComments );
      area.setEditable(false);
      panel_comments.add(area, BorderLayout.NORTH);
    }

    JPanel panelAll = new JPanel(new BorderLayout());
    panelAll.setBorder( BorderFactory.createCompoundBorder(
        BorderFactory.createEmptyBorder( 8, 11, 11, 11 ), csStandard.INNER_PARAM_BORDER) );
    panelAll.add( panelIcon, BorderLayout.NORTH );
    panelAll.add( panel, BorderLayout.CENTER );
    panelAll.add( panel_comments, BorderLayout.SOUTH );
    
    panel.setBackground( Color.white );
    panelIcon.setBackground( Color.white );
    panelAll.setBackground( Color.white );

    getContentPane().add(panelAll,BorderLayout.CENTER);
    pack();
    setLocationRelativeTo( myParentFrame );
  }
  public void setTitle( String title ) {
    setTitle( title );
  }
  public void setLogo( ImageIcon logo ) {
    myLogo = logo;
  }
  public void setVersionString( String versionString ) {
    myVersion = versionString;
  }
  public void setAuthor( String author ) {
	  myAuthor = author;
	}
  public void setAdditionalComments( String comments ) {
	  myComments = comments;
	}
  public void setProductName( String productName ) {
    myProductName = productName;
  }
  public void setDate( String dateString ) {
    myDate = dateString;
  }
  public void setCompany( String company ) {
    myCompany = company;
  }
  public void setContact( String contact ) {
    myContact = contact;
  }

}


