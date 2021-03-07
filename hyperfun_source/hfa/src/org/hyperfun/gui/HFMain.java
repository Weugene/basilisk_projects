/*
 * Copyright (c) 2002-2003 Mio Hiraga, Yuichiro Goto
 * All rights reserved.
 *
 * This program can be redistributed and/or modified under the terms
 * of the CGPL, The Common Good Public License as published by and at CGPL.org
 * (http://CGPL.org).  It is released under version 1.0 Beta of the License
 * until the 1.0 version is released after which either version 1.0 of the
 * License, or (at your option) any later version can be applied.
 *
 * THIS WORK, OR SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED (See the
 * CGPL, The Common Good Public License for more information.)
 *
 * You should have received a copy of the CGPL along with this program;
 * if not, see -  http://CGPL.org to get a copy of the License.
 *
 * $Id: HFMain.java,v 1.7 2006/11/08 09:23:49 y7goto Exp $
 */

/*
 * modified PA 2 11 2004 added strings for convolution primitives
 * modified PA 2 12 2004 removed description of hfConvLineR
 * modified the bounding box for the object
 */

package org.hyperfun.gui;

import java.util.ArrayList;
import java.util.Hashtable;
import java.io.*;

import java.awt.*;
import java.awt.event.*;
import java.applet.*;
import javax.swing.*;
import javax.swing.text.*;
import javax.swing.border.*;

import javax.media.j3d.*;
import javax.vecmath.*;
import com.sun.j3d.utils.universe.SimpleUniverse;
import com.sun.j3d.utils.behaviors.mouse.MouseRotate;
import com.sun.j3d.utils.behaviors.mouse.MouseTranslate;
import com.sun.j3d.utils.behaviors.mouse.MouseZoom;
import com.sun.j3d.utils.geometry.Text2D;

import org.hyperfun.EmpiricalHyperFun.Calculator;
import org.hyperfun.polygonizer.Polygonizer;

/**
 * GUI part of the HyperFun polygonal viewer.  Since this class has the main
 * method, it can run as both the stand-alone Java application and Applet.  In
 * case running as the Applet, the jar file has to be signed using jarsigner,
 * because the applet uses a custom class loader, and it accesses local files
 * when the resulting polygonal model is saved.
 *
 * @author Mio Hiraga
 * @author Yuichiro Goto
 * @version 0.2
 */
public class HFMain extends JApplet {
    private static final int X = 0;
    private static final int Y = 1;
    private static final int Z = 2;

    // For drawing Model by Java 3D
    Shape3D object = new Shape3D();

    // For drawing BoundingBox by Java 3D
    Shape3D boundingBox = new Shape3D();

    // For drawing XYZ-axis by Java 3D
    Shape3D axis = new Shape3D();

    // For setting geometry datas on Model by Java 3D
    GeometryArray triangles = null;

    // For setting geometry datas on BoundingBox by Java 3D
    GeometryArray bdbox = null;

    // For setting geometry datas on BoundingBox by Java 3D
    GeometryArray axis_geom = null;

    Background background = new Background();

    // Vertices datas of Model
    Point3d[] vertices;

    // Normals datas of Model
    Vector3f[] normals;

    // Indecices datas of Model
    int[] indices;

    // Default color of Model
    Color ObjectColor = new Color(0.0f, 0.8f, 0.2f);

    // For showing on the applet
    String defaultLowerBounds = "-5.0";

    // For showing on the applet
    String defaultUpperBounds = "5.0";

    // For showing on the applet
    String defaultGrid = "30";

    // For easy test
    //String defaultGrid = "10";

    // For showing output text on the save window
    String outputtxt = "";

    boolean isJapanese = false;
    boolean isEnglish = true;
    boolean toJapanese = false;
    boolean toEnglish = true;
    boolean isPolygonizedButtonPushed = false;
    boolean isHfSaveButtonPushed = false;
    boolean isWrlSaveButtonPushed = false;
    boolean isStlSaveButtonPushed = false;
    boolean isStandalone = false;

    Polygonizer polygonizer = new Polygonizer(
	new double[] {-5.0, -5.0, -5.0},
	// Bounding Box Value
	new double[] {5.0, 5.0, 5.0},
	// Grid Size, IsoValue
	new int[] {30, 30, 30}, 0.0,
	null);

    //*****GUI by Swing of Java2**************************************************

    // Base of the panel
    JPanel panelBase = new JPanel();

    // For splitting half left and right
    JSplitPane jSplitPane1 = new JSplitPane();

    //***right panel***
    JPanel panelRight = new JPanel();

    // Textarea for hfcode
    JTextArea jTextArea1 = new JTextArea();

    // Making scroll bar for textarea
    JScrollPane jScrollPane1 = new JScrollPane();

    // Tab menu
    JTabbedPane tabbedPane = new JTabbedPane();

    // Tab menu for BoundingBox
    JPanel panel1 = new JPanel();

    // Tab menu for Grid
    JPanel panel2 = new JPanel();

    // Tab menu for Color
    JPanel panel3 = new JPanel();

    // Tab menu for Library
    JPanel panel4 = new JPanel();

    // Tab menu for Save
    JPanel panel5 = new JPanel();

    // Panel for Button
    JPanel panelB = new JPanel();

    // Panel for BoundingBox menu
    JPanel panel_bdbox = new JPanel();

    // Label for BoundingBox menu
    JLabel lb_xmin, lb_ymin, lb_zmin, lb_xmax, lb_ymax, lb_zmax;

    // Input textarea for BoundingBox menu
    JTextField txt_xmin, txt_ymin, txt_zmin, txt_xmax, txt_ymax, txt_zmax;

    // Reset button for Bounding Box values
    JButton resetBdBoxButton = new JButton();

    // Panel for Grid menu
    JPanel panel_grid = new JPanel();

    // Label for Grid menu
    JLabel lb_gridx, lb_gridy, lb_gridz;

    // Label for Grid menu
    JTextField txt_gridx, txt_gridy, txt_gridz;

    // Reset button for Grid values
    JButton resetGridButton = new JButton();

    // Panel for Color menu
    JPanel panelColor = new JPanel();

    // Select Button for Color
    JButton selectButton = new JButton();

    // Label for Color menu
    static JLabel labelColor = new JLabel();

    // Selecting Model color
    static JColorChooser colorchooser = new JColorChooser();

    // Japanese Button
    JButton jpnButton = new JButton();

    // PolygonizeButton
    JButton polygonizeButton = new JButton();

    // SaveButton as HyperFun
    JButton hfSaveButton = new JButton();

    // SaveButton as VRML
    JButton wrlSaveButton = new JButton();

    // SaveButton as STL
    JButton stlSaveButton = new JButton();

    JButton povraySaveButton = new JButton();

    // File choose window
    JFileChooser fileChooser = new JFileChooser();

    //***left panel***
    JPanel panelLeft = new JPanel();

    // For splitting top and bottom on the left panel
    JSplitPane jSplitPane2 = new JSplitPane();

    // Textarea for error messages
    JTextArea errorTextArea = new JTextArea();

    // Making scroll bar for error textarea
    JScrollPane jScrollPane2 = new JScrollPane();

    // Pop up menu on the textarea hor hfcode
    JPopupMenu popup = new JPopupMenu();

    // Pop up window for showing file
    JFrame fileViewer = new JFrame();

    JButton hfOpenButton = new JButton();

    Hashtable actions = null;

    String defaultHFCode =
	"my_model(x[3], a[1])\n" +
	"{\n" +
	" array center[3];\n" +
	" center = [0,2,0];\n"   +
	" torus = hfTorusZ(x,center,2,0.5);\n" +
	" xx = x[1]; yy = x[2]; zz = x[3];\n"+
	" piriform = -( (yy/5)^4+(yy/5)^3\n" +
	" 		+(zz/7)^2+(xx/7)^2 );\n" +
	"\n" +
	" my_model = torus | piriform ;\n" +
	"}\n";

    public String getParameter(String key, String def) {
	return isStandalone ? System.getProperty(key, def) : (getParameter(key) != null ? getParameter(key) : def);
    }

    public HFMain() {
    }

    public void init() {
	try {
	    jbInit();
	} catch(Exception e) {
	    e.printStackTrace();
	}
    }

    private void jbInit() throws Exception {
	this.setSize(new Dimension(640, 480));

	GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
	Canvas3D canvas = new Canvas3D(config);
	SimpleUniverse universe = new SimpleUniverse(canvas);
	BranchGroup scene = createSceneGraph();
	universe.addBranchGraph(scene);
        universe.getViewer().getView().setBackClipDistance(100.0);

	this.getContentPane().add(panelBase, BorderLayout.CENTER);

	// Splitting left and right
	jSplitPane1.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
	jSplitPane1.setDividerLocation(320);
	panelBase.setLayout(new BorderLayout());
	panelBase.add(jSplitPane1, BorderLayout.CENTER);

	//***********RIGHT Side*******************************************************

	panelRight.setSize(320, 480);
	panelRight.setLayout(new BorderLayout(0, 8));
	panelRight.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
	// Setting on the right side
	jSplitPane1.setRightComponent(panelRight);

	//***HyperFunCode
	jScrollPane1.setPreferredSize(new Dimension(280, 220));
	jScrollPane1.getViewport().add(jTextArea1, null);
	jTextArea1.setText(defaultHFCode);
	jTextArea1.setCaretPosition(0);
	jTextArea1.setMargin(new Insets(5, 5, 5, 5));
	panelRight.add(jScrollPane1, BorderLayout.NORTH);

	//***Tab menu***********
	// Icon on the tab menu
	ImageIcon icon = null;
	tabbedPane.addTab("Model Space", icon, panel1, "You can change your model space.");
	tabbedPane.setSelectedIndex(0);
	panel1.setLayout(new BorderLayout(10, 10));
	panel1.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));

	//**Bounding Box menu
	TitledBorder title;
	title = BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Enter the bounding box values.");
	panel_bdbox.setBorder(title);
	panel_bdbox.setLayout(new GridLayout(2, 6));
	panel1.add(panel_bdbox, BorderLayout.NORTH);

	lb_xmin = new JLabel("Xmin: ");
	lb_xmin.setHorizontalAlignment(JLabel.RIGHT);
	panel_bdbox.add(lb_xmin);
	txt_xmin = new JTextField(defaultLowerBounds, 4);
	panel_bdbox.add(txt_xmin);

	lb_ymin = new JLabel("Ymin: ");
	lb_ymin.setHorizontalAlignment(JLabel.RIGHT);
	panel_bdbox.add(lb_ymin);
	txt_ymin = new JTextField(defaultLowerBounds, 4);
	panel_bdbox.add(txt_ymin);

	lb_zmin = new JLabel("Zmin: ");
	lb_zmin.setHorizontalAlignment(JLabel.RIGHT);
	panel_bdbox.add(lb_zmin);
	txt_zmin = new JTextField(defaultLowerBounds, 4);
	panel_bdbox.add(txt_zmin);

	lb_xmax = new JLabel("Xmax: ");
	lb_xmax.setHorizontalAlignment(JLabel.RIGHT);
	panel_bdbox.add(lb_xmax);
	txt_xmax = new JTextField(defaultUpperBounds, 4);
	panel_bdbox.add(txt_xmax);

	lb_ymax = new JLabel("Ymax: ");
	lb_ymax.setHorizontalAlignment(JLabel.RIGHT);
	panel_bdbox.add(lb_ymax);
	txt_ymax = new JTextField(defaultUpperBounds, 4);
	panel_bdbox.add(txt_ymax);

	lb_zmax = new JLabel("Zmax: ");
	lb_zmax.setHorizontalAlignment(JLabel.RIGHT);
	panel_bdbox.add(lb_zmax);
	txt_zmax = new JTextField(defaultUpperBounds, 4);
	panel_bdbox.add(txt_zmax);

	JPanel buttonPane = new JPanel();
	buttonPane.setLayout(new BoxLayout(buttonPane, BoxLayout.X_AXIS));
	buttonPane.add(Box.createHorizontalGlue());
	buttonPane.add(resetBdBoxButton);

	resetBdBoxButton.setText("reset");
	resetBdBoxButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {resetBdBoxButton_actionPerformed(e);}
	    }
	);
	panel1.add(buttonPane, BorderLayout.CENTER);

	//***Grid menu
	tabbedPane.addTab("Grid", icon, panel2, "You can change the grid size.");
	panel2.setLayout(new BorderLayout(10, 10));
	panel2.setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

	//panel_grid.setPreferredSize(new Dimension(180, 100));
	panel_grid.setLayout(new GridLayout(3, 2));
	title = BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Enter the number of grid lines.");
	panel_grid.setBorder(title);

	lb_gridx = new JLabel("Grid X :  ");
	txt_gridx = new JTextField(defaultGrid, 4);
	lb_gridy = new JLabel("Grid Y :  ");
	txt_gridy = new JTextField(defaultGrid, 4);
	lb_gridz = new JLabel("Grid Z :  ");
	txt_gridz = new JTextField(defaultGrid, 4);
	lb_gridx.setHorizontalAlignment(JLabel.RIGHT);
	panel_grid.add(lb_gridx);
	panel_grid.add(txt_gridx);
	lb_gridy.setHorizontalAlignment(JLabel.RIGHT);
	panel_grid.add(lb_gridy);
	panel_grid.add(txt_gridy);
	lb_gridz.setHorizontalAlignment(JLabel.RIGHT);
	panel_grid.add(lb_gridz);
	panel_grid.add(txt_gridz);
	panel2.add(panel_grid,  BorderLayout.NORTH);

	resetGridButton.setText("reset");
	resetGridButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {resetGridButton_actionPerformed(e);}
	    }
	);

	JPanel buttonPane2 = new JPanel();
	buttonPane2.setLayout(new BoxLayout(buttonPane2, BoxLayout.X_AXIS));
	buttonPane2.add(Box.createHorizontalGlue());
	buttonPane2.add(resetGridButton);
	panel2.add(buttonPane2,  BorderLayout.CENTER);

	tabbedPane.addTab("Color", icon, panel3, "You can change your object color.");
	panel3.setLayout(new BorderLayout(10, 10));
	panel3.setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));

	title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Select the color of your object.");
	panelColor.setBorder(title);

	labelColor.setPreferredSize(new Dimension(60, 60));
	labelColor.setOpaque(true);
	// Default object color
	labelColor.setBackground(ObjectColor);
	labelColor.setBorder(BorderFactory.createLineBorder(Color.black));
	panelColor.add(labelColor);
	panel3.add(panelColor, BorderLayout.NORTH);

	selectButton.setText("select");
	selectButton.addActionListener(
	    new ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    Color color = colorchooser.showDialog(labelColor, "Choose your color !", ObjectColor);
		    //System.out.println("R: "+color.getRed() + "G: "+color.getGreen() + "B: "+color.getBlue());
		    if (color != null) {
			labelColor.setBackground( color );
			ObjectColor = color;
		    }
		}
	    }
	);

	panelColor.add(selectButton);

	//***Library
	tabbedPane.addTab("Library", icon, panel4, "You can reference HF library.");

	String lib =
	    "<html><body><p><b>Primitives:</b><br>" +
	    " - hfSphere(x, center, R);<br>" +
	    " - hfEllipsoid(x, center, a, b, c);<br>" +
	    " - hfCylinderX(x, center,R);<br>" +
	    " - hfCylinderY(x, center,R);<br>" +
	    " - hfCylinderZ(x, center,R);<br>" +
	    " - hfEllCylX(x, center, a, b);<br>" +
	    " - hfEllCylY(x, center, a, b);<br>" +
	    " - hfEllCylZ(x, center, a, b);<br>" +
	    " - hfTorusX(x, center, R, r0);<br>" +
	    " - hfTorusY(x, center, R, r0);<br>" +
	    " - hfTorusZ(x, center, R, r0);<br>" +
	    " - hfBlock(x, vertex, dx, dy, dz);<br>" +
	    " - hfBlobby(x,x0,y0,z0,a,b,T);<br>" +
	    " - hfMetaBall(x,x0,y0,z0,b,d,T);<br>" +
	    " - hfSoft(x,x0,y0,z0,d,T);<br>" +
		/*New primitives BEGIN*/
	    " - hfConvArc(x, center, radius, theta, axis, angle, S, T);<br>" +
	    " - hfConvCurve(x, vect, S, T);<br>" +
	    " - hfConvLine(x, begin, end, S, T);<br>" +
	    " - hfConvMesh(x, vect, tri, S, T);<br>" +
	    " - hfConvTriangle(x, vect, S, T);<br>" +
	    " - hfConvPoint(x, vect, S, T);<br>" +
	    " - hfConeX(x, center, R);<br>" +
	    " - hfConeY(x, center, R);<br>" +
	    " - hfConeZ(x, center, R);<br>" +
	    " - hfEllConeX(x, center, a, b);<br>" +
	    " - hfEllConeY(x, center, a, b);<br>" +
	    " - hfEllConeZ(x, center, a, b);<br>" +
	    " - hfNoiseG(x, amp, freq, phase);<br>" +
	    " - hfSuperEll(x, center, a, b, c, s1, s2);</p>" +
	    /*New primitives END*/
	    "<p><br><b>Operations:</b><br>" +
	    " - hfScale3D(xt, sx, sy, sz);<br>" +
	    " - hfShift3D(xt, dx, dy, dz);<br>" +
	    " - hfRotate3DX(xt, theta);<br>" +
	    " - hfRotate3DY(xt, theta);<br>" +
	    " - hfRotate3DZ(xt, theta);<br>" +
	    " - hfBlendUni(f1,f2,a0,a1,a2);<br>" +
	    " - hfBlendInt(f1,f2,a0,a1,a2);<br>" +
	    " - hfTwistX(xt, x1, x2, theta1, theta2);<br>" +
	    " - hfTwistY(xt, y1, y2, theta1, theta2);<br>" +
	    " - hfTwistZ(xt, z1, z2, theta1, theta2);<br>" +
	    " - hfStretch3D(xt,x0,sx,sy,sz);<br>" +
	    " - hfTaperX(xt, x1, x2, s1, s2);<br>" +
	    " - hfTaperY(xt, y1, y2, s1, s2);<br>" +
	    " - hfTaperZ(xt, z1, z2, s1, s2);<br>" +
	    " - hfFMapBlob(x, fobj, x0, y0, z0, fobj0, sigma);<br>" +
	    " - hfSpaceMapCubic(x, source, target, b);<br>" +
	    " - hfSpaceMapExp(x, source, target, b);</p></body></html>";
	JLabel libtxt = new JLabel(lib);
	JScrollPane jScrollPane3 = new JScrollPane();
	jScrollPane3.setPreferredSize(new Dimension(260, 140));
	jScrollPane3.getViewport().add(libtxt, null);
	panel4.add(jScrollPane3);

	//***Save
	tabbedPane.addTab("Save", icon, panel5, "You can save the file.");
	panel5.setLayout(new BorderLayout());

	//***save as HF
	// Space(top, left, bottom, right)
	panel5.setBorder(BorderFactory.createEmptyBorder(20, 60, 20, 60));
	hfSaveButton.setText("Save as HF");
	hfSaveButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    hfSaveButton_actionPerformed(e);
		}
	    }
	);
	JPanel buttonPane5 = new JPanel();
	buttonPane5.setLayout(new GridLayout(5, 1, 10, 10));
	buttonPane5.add(hfSaveButton);

	// Save as hf file
	hfOpenButton.setText("Open hf file");
	hfOpenButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    openButton_actionPerformed(e);
		}
	    }
	);
	buttonPane5.add(hfOpenButton);

	//***save as VRML
	wrlSaveButton.setText("Save as VRML");
	wrlSaveButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    wrlSaveButton_actionPerformed(e);
		}
	    }
	);
	buttonPane5.add(wrlSaveButton);

	//***save as STL
	stlSaveButton.setText("Save as STL");
	stlSaveButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    stlSaveButton_actionPerformed(e);
		}
	    }
	);
	buttonPane5.add(stlSaveButton);


	// Save mesh povray
	povraySaveButton.setText("Save as povray mesh");
	povraySaveButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e){
		    povraySaveButton_actionPerformed(e);
		}
	    }
	);

	buttonPane5.add(povraySaveButton);

	panel5.add(buttonPane5,  BorderLayout.CENTER);

	// Set tab menu on the right panel
	panelRight.add(tabbedPane, BorderLayout.CENTER);

	panelB.setLayout(new BorderLayout(10, 10));
	//***Japanese Button
	jpnButton.setText("Japanese");
	jpnButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {
		    if (isJapanese == true) {
			toEnglish = true;
		    }
		    if (isEnglish == true) {
			toJapanese = true;
		    }
		    jpnButton_actionPerformed(e);
		}
	    }
	);
	panelB.add(jpnButton,  BorderLayout.WEST);

	//***Polygonize Button
	polygonizeButton.setText(" POLYGONIZE ");
	polygonizeButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {polygonizeButton_actionPerformed(e);}
	    }
	);
	panelB.add(polygonizeButton,  BorderLayout.CENTER);

	panelRight.add(panelB, BorderLayout.SOUTH);

	//***********LEFT Side*******************************************************

	panelLeft.setPreferredSize(new Dimension(320, 480));
	panelLeft.setLayout(new BorderLayout());
	jSplitPane1.setLeftComponent(panelLeft);

	// Splitting top and bottom
	jSplitPane2.setOrientation(JSplitPane.VERTICAL_SPLIT);
	jSplitPane2.setDividerLocation(320);
	panelLeft.add(jSplitPane2, BorderLayout.CENTER);

	jSplitPane2.setTopComponent(canvas);

	errorTextArea.setEditable(false);
	errorTextArea.setCaretPosition(0);
	jScrollPane2.getViewport().add(errorTextArea, null);
	jSplitPane2.setBottomComponent(jScrollPane2);

	/***************************************************************************/

	// Create a popup menu for editing
	createActionTable(jTextArea1);

	JMenuItem item = new JMenuItem();
	item.setAction(getActionByName(DefaultEditorKit.cutAction));
	popup.add(item);
	item = new JMenuItem();
	item.setAction(getActionByName(DefaultEditorKit.copyAction));
	popup.add(item);
	item = new JMenuItem();
	item.setAction(getActionByName(DefaultEditorKit.pasteAction));
	popup.add(item);
	popup.add(new JSeparator());
	item = new JMenuItem();
	item.setAction(getActionByName(DefaultEditorKit.selectAllAction));
	popup.add(item);

	MouseListener popupListener = new PopupListener();
	jTextArea1.addMouseListener(popupListener);
    }

    //*****Create Scene by Java 3D***************************************************
    private BranchGroup createSceneGraph() {
	BranchGroup root = new BranchGroup();
	//root.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);

	BoundingBox bounds = new BoundingBox(new Point3d(-200.0, -200.0, -200.0), new Point3d(200.0, 200.0, 200.0));

	//Background background = new Background();
	background.setApplicationBounds(bounds);
	background.setColor(new Color3f(0.5f, 0.5f, 0.5f));
	root.addChild(background);

	Transform3D t3d = new Transform3D();
	t3d.rotX(-Math.PI/2.0);
	t3d.setTranslation(new Vector3d(0.0, 0.0, -15.0));
	TransformGroup trans = new TransformGroup(t3d);

	trans.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
	trans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);

	MouseRotate rotate = new MouseRotate(trans);
	rotate.setSchedulingBounds(bounds);
	root.addChild(rotate);
	MouseTranslate translate = new MouseTranslate(trans);
	translate.setSchedulingBounds(bounds);
	root.addChild(translate);
	MouseZoom zoom = new MouseZoom(trans);
	zoom.setSchedulingBounds(bounds);
	root.addChild(zoom);

	DirectionalLight dlight = new DirectionalLight(
	    new Color3f(1.0f, 1.0f, 1.0f),
	    new Vector3f(0.0f, 0.0f, -1.0f)
	);
	dlight.setInfluencingBounds(bounds);
	root.addChild(dlight);

	object.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
	object.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
	boundingBox.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
	axis.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);

	trans.addChild(object);
	trans.addChild(boundingBox);
	trans.addChild(axis);
	root.addChild(trans);

	return root;
    }

    private void setMaterial(Color color) {
	Material mat = new Material(
	    // Ambient color
	    new Color3f(0.2f, 0.2f, 0.2f),
	    // Diffuse color
	    new Color3f(0.0f, 0.0f, 0.0f),
	    // Specular color
	    new Color3f(color),
	    // Emissive color
	    new Color3f(0.0f, 0.0f, 0.0f),
	    // Shininess
	    64.0f);

	PolygonAttributes poly = new PolygonAttributes(
	    PolygonAttributes.POLYGON_FILL,
	    PolygonAttributes.CULL_NONE,
	    0.0f);
	Appearance app = new Appearance();
	app.setMaterial(mat);
	app.setPolygonAttributes(poly);
	object.setAppearance(app);
    }

    private GeometryArray createAxis(double[] lowerBounds, double[] upperBounds) {
	// axis
	Point3d[] axis_vertices = new Point3d[6];

	axis_vertices[0] = new Point3d(0.0, 0.0, 0.0);
	axis_vertices[1] = new Point3d(upperBounds[0], 0.0, 0.0);
	axis_vertices[2] = new Point3d(0.0, 0.0, 0.0);
	axis_vertices[3] = new Point3d(0.0, upperBounds[1], 0.0);
	axis_vertices[4] = new Point3d(0.0, 0.0, 0.0);
	axis_vertices[5] = new Point3d(0.0, 0.0, upperBounds[2]);

	LineArray axis_geometry = new LineArray(axis_vertices.length, GeometryArray.COORDINATES | GeometryArray.COLOR_3);
	axis_geometry.setCoordinates(0, axis_vertices);
	axis_geometry.setColor(0, new Color3f(1.0f, 0.0f, 0.0f));
	axis_geometry.setColor(1, new Color3f(1.0f, 0.0f, 0.0f));
	axis_geometry.setColor(2, new Color3f(0.0f, 1.0f, 0.0f));
	axis_geometry.setColor(3, new Color3f(0.0f, 1.0f, 0.0f));
	axis_geometry.setColor(4, new Color3f(0.0f, 0.0f, 1.0f));
	axis_geometry.setColor(5, new Color3f(0.0f, 0.0f, 1.0f));

	return axis_geometry;
    }

    private GeometryArray createBoundingBox(double[] lowerBounds, double[] upperBounds) {
	Point3d[] vertices = {
	    new Point3d(lowerBounds[0], lowerBounds[1], lowerBounds[2]),
	    new Point3d(lowerBounds[0], upperBounds[1], lowerBounds[2]),
	    new Point3d(upperBounds[0], upperBounds[1], lowerBounds[2]),
	    new Point3d(upperBounds[0], lowerBounds[1], lowerBounds[2]),
	    new Point3d(lowerBounds[0], lowerBounds[1], upperBounds[2]),
	    new Point3d(lowerBounds[0], upperBounds[1], upperBounds[2]),
	    new Point3d(upperBounds[0], upperBounds[1], upperBounds[2]),
	    new Point3d(upperBounds[0], lowerBounds[1], upperBounds[2])
	};

	int[] indices = {
	    0, 1,
	    1, 2,
	    2, 3,
	    3, 0,
	    4, 5,
	    5, 6,
	    6, 7,
	    7, 4,
	    0, 4,
	    1, 5,
	    2, 6,
	    3, 7
	};

	IndexedLineArray geom = new IndexedLineArray(
	    vertices.length,
	    GeometryArray.COORDINATES,
	    indices.length);
	geom.setCoordinates(0, vertices);
	geom.setCoordinateIndices(0, indices);

	return geom;
    }

    private GeometryArray getTriangles() {
	/*
	Point3d[] vertices;
	Vector3f[] normals;
	Color3f[] colors;
	int[] indices;
	*/
	ArrayList polygon_vertices = new ArrayList();
	ArrayList polygon_normals = new ArrayList();
	ArrayList polygon_indices = new ArrayList();
	int i;

	try {
	    polygonizer.polygonize(polygon_vertices,
				   polygon_normals,
				   polygon_indices);
	} catch (Exception e) {
	    System.out.println(e);
	    return null;
	}

	if (polygon_vertices.size() == 0) {
	    return null;
	}

	vertices = new Point3d[polygon_vertices.size()];
	normals = new Vector3f[polygon_normals.size()];
	indices = new int[polygon_indices.size()];
	//colors = new Color3f[polygon_vertices.size()];

	java.util.Iterator vitr;
	java.util.Iterator nitr;

	for (vitr = polygon_vertices.iterator(), i = 0, nitr = polygon_normals.iterator();
	     vitr.hasNext() && nitr.hasNext(); i++) {
	    vertices[i] = new Point3d((double[])vitr.next());
	    normals[i] = new Vector3f((float[])nitr.next());
	    //colors[i] = new Color3f(0.0f, 0.5f, 0.0f);
	}

	java.util.Iterator iitr;
	for (iitr = polygon_indices.iterator(), i = 0; iitr.hasNext(); i++) {
	    indices[i] = ((Integer)iitr.next()).intValue();
	    //System.out.println("index: " + indices[i]);
	}

	/*////// Print debug information///////////
	System.out.println( "length:" + vertices.length );
	for (int ii = 0; ii < vertices.length; ii++) {
	    int nn = 0;
	    System.out.println("vertices[" + ii + "] = " + vertices[ii]);
	    for (int jj = ii+1; jj < vertices.length; jj++) {
		if (vertices[ii].distance( vertices[jj] ) < 0.000001) {
		    nn++;
		    System.out.println("here:" + nn);
		}
	    }
	}
	for (int iii = 0; iii < normals.length; iii++) {
	    System.out.println("normals[" + iii + "] = " + normals[iii]);
	}
	for (int iii = 0; iii < indices.length; iii++) {
	    System.out.println("indices[" + iii + "] = " + indices[iii]);
	}
	//////////////////////////////////////*/

	IndexedTriangleArray geom = new IndexedTriangleArray(
	    vertices.length,
	    GeometryArray.COORDINATES | GeometryArray.NORMALS,
	    indices.length);
	geom.setCoordinates(0, vertices);
	geom.setCoordinateIndices(0, indices);
	geom.setNormals(0, normals);
	geom.setNormalIndices(0, indices);
	//geom.setColors(0, colors);
	//geom.setColorIndices(0, indices);

	return geom;
    }

    public void start() {
    }

    public void stop() {
    }

    public void destroy() {
    }

    public String getAppletInfo() {
	return "Applet Information";
    }

    public String[][] getParameterInfo() {
	return null;
    }

    public static void main(String[] args) {
	HFMain applet = new HFMain();
	applet.isStandalone = true;

	JFrame frame = new JFrame();
	frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	frame.setTitle("Hyper FUNFUN");
	frame.getContentPane().add(applet, BorderLayout.CENTER);

	applet.init();
	applet.start();
	frame.setSize(640, 480);

	Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
	frame.setLocation((d.width - frame.getSize().width)/2, (d.height - frame.getSize().height)/2);
	frame.setVisible(true);
    }

    static {
	try {
	  UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
	  //UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
	} catch (Exception e) {
	}
    }

    private void createActionTable(JTextComponent textComponent) {
	actions = new Hashtable();
	Action[] actionsArray = textComponent.getActions();

	for (int i = 0; i < actionsArray.length; i++) {
	    Action a = actionsArray[i];
	    actions.put(a.getValue(Action.NAME), a);
	}
    }

    private Action getActionByName(String name) {
	return (Action)(actions.get(name));
    }

    // To trigger popup menu
    class PopupListener extends MouseAdapter {
	public void mousePressed(MouseEvent e) {
	    showPopup(e);
	}
	public void mouseReleased(MouseEvent e) {
	    showPopup(e);
	}
	private void showPopup(MouseEvent e) {
	    if (e.isPopupTrigger()) {
		popup.show(e.getComponent(), e.getX(), e.getY());
	    }
	}
    }

    void polygonizeButton_actionPerformed(ActionEvent event) {
	double[] upperBounds = new double[3];
	double[] lowerBounds = new double[3];
	int[] div = new int[3];

	String buf = jTextArea1.getText();

	// buf to lower case
	// hf specification:
	// code is case insensitive
	buf = buf.toLowerCase();

	if (buf.length() != 0) {
	    try {
		// Clear previous error messages on errorTextArea
		errorTextArea.setText(null);
		isPolygonizedButtonPushed = true;

		lowerBounds[X] = Double.parseDouble(txt_xmin.getText());
		lowerBounds[Y] = Double.parseDouble(txt_ymin.getText());
		lowerBounds[Z] = Double.parseDouble(txt_zmin.getText());
		upperBounds[X] = Double.parseDouble(txt_xmax.getText());
		upperBounds[Y] = Double.parseDouble(txt_ymax.getText());
		upperBounds[Z] = Double.parseDouble(txt_zmax.getText());
		polygonizer.setBoundingBox(lowerBounds, upperBounds);
		polygonizer.getBoundingBox(lowerBounds, upperBounds);
		bdbox = createBoundingBox(lowerBounds, upperBounds);
		// BoundingBox drawing
		boundingBox.setGeometry(bdbox);

		axis_geom = createAxis(lowerBounds, upperBounds);
		// XYZ-axis drawing
		axis.setGeometry(axis_geom);

		div[X] = Integer.parseInt(txt_gridx.getText());
		div[Y] = Integer.parseInt(txt_gridy.getText());
		div[Z] = Integer.parseInt(txt_gridz.getText());
		// Set Grid Size
		polygonizer.setDivisions(div[X], div[Y], div[Z]);

		Calculator calc = new Calculator(buf);
		polygonizer.setFunction(calc);
		triangles = getTriangles();

		if (triangles != null) {
		    setMaterial(ObjectColor);
		    // Object drawing
		    object.setGeometry(triangles);
		}
	    } catch (Exception exception) {
		errorTextArea.setText("Parse error\n");
		errorTextArea.append(exception.getMessage());
		//System.out.println(exception.getClass().getName() + ": " + exception.getMessage());
	    }
	}
    }

    void hfSaveButton_actionPerformed(ActionEvent event) {
	isHfSaveButtonPushed = true;
	outputtxt = jTextArea1.getText();
	showOutPutWindow();
    }

    void wrlSaveButton_actionPerformed(ActionEvent event) {
	isWrlSaveButtonPushed = true;

	if (isPolygonizedButtonPushed == true) {
	    String buf = jTextArea1.getText();
	    VrmlOut vrmlout = new VrmlOut(buf , ObjectColor , vertices, normals, indices);
	    outputtxt =  vrmlout.vrmlout();
	    showOutPutWindow();
	} else {
	    errorTextArea.setText("Please save the file after POLYGONIZE ! \n");
	}
    }

    void stlSaveButton_actionPerformed(ActionEvent event) {
	isStlSaveButtonPushed = true;

	if (isPolygonizedButtonPushed == true) {
	    StlOut stlout = new StlOut(vertices, indices);

	    ExampleFileFilter stlFilter = new ExampleFileFilter("stl", "stl mesh files");

	    File stl = new File("files/output.stl");

	    fileChooser.addChoosableFileFilter(stlFilter);
	    fileChooser.setSelectedFile(stl);

	    fileChooser.setAcceptAllFileFilterUsed(false);
	    int returnVal = fileChooser.showSaveDialog(this);

	    try{
		if(returnVal == JFileChooser.APPROVE_OPTION) {
		    File file = fileChooser.getSelectedFile();
		    stlout.write(file);
		}

	    } catch (Exception ex){
		ex.printStackTrace();
	    }
	} else {
	    errorTextArea.setText("Please save the file after POLYGONIZE ! \n");
       	}
    }

    void povraySaveButton_actionPerformed(ActionEvent event){
	if (isPolygonizedButtonPushed == true) {
	    double[] l = new double[3];
	    double[] u = new double[3];
	    polygonizer.getBoundingBox(l,u);
	    PovrayOut out = new PovrayOut(vertices, normals, indices, l, u);

	    ExampleFileFilter povrayFilter = new ExampleFileFilter("pov", "povray script files");

	    File povray = new File("files/output.pov");

	    fileChooser.addChoosableFileFilter(povrayFilter);
	    fileChooser.setSelectedFile(povray);

	    fileChooser.setAcceptAllFileFilterUsed(false);
	    int returnVal = fileChooser.showSaveDialog(this);

	    try{
		if (returnVal == JFileChooser.APPROVE_OPTION) {
		    File file = fileChooser.getSelectedFile();
		    out.write(file);
		}
	    } catch (Exception ex){
		ex.printStackTrace();
	    }
	} else {
	    errorTextArea.setText("Please save the file after POLYGONIZE !\n");
	}
    }

    void resetBdBoxButton_actionPerformed(ActionEvent event){
	txt_xmin.setText(defaultLowerBounds);
	txt_ymin.setText(defaultLowerBounds);
	txt_zmin.setText(defaultLowerBounds);
	txt_xmax.setText(defaultUpperBounds);
	txt_ymax.setText(defaultUpperBounds);
	txt_zmax.setText(defaultUpperBounds);
    }

    void resetGridButton_actionPerformed(ActionEvent event){
	txt_gridx.setText(defaultGrid);
	txt_gridy.setText(defaultGrid);
	txt_gridz.setText(defaultGrid);
    }

    void showOutPutWindow() {
	//***show the output window for saving a file

	JFrame outputWindow = new JFrame("output file");
	JPanel buttonPanel = new JPanel();
	JLabel label = new JLabel();
	JButton saveButton = new JButton("Save");
	JScrollPane outputjScrollPane = new JScrollPane();
	JTextArea outputTextArea = new JTextArea();

	if (isHfSaveButtonPushed == true) {
	    label.setText("HyperFun Format   ");
	}
	if (isWrlSaveButtonPushed == true) {
	    label.setText("VRML Format   ");
	}
	if (isStlSaveButtonPushed == true) {
	    label.setText("StL ASCII Format   ");
	}

	buttonPanel.add(label);
	buttonPanel.add(saveButton);
	outputWindow.getContentPane().add(buttonPanel, BorderLayout.NORTH);
	outputjScrollPane.setPreferredSize(new Dimension(300,400));
	outputjScrollPane.getViewport().add(outputTextArea, null);
	outputTextArea.setText(outputtxt);
	outputTextArea.setCaretPosition(0);
	outputTextArea.setMargin(new Insets(5, 5, 5, 5));
	outputWindow.getContentPane().add(outputjScrollPane, BorderLayout.CENTER);

	saveButton.addActionListener(
	    new java.awt.event.ActionListener() {
		public void actionPerformed(ActionEvent e) {saveButton_actionPerformed(e);}
	    }
	);
	outputWindow.setBounds( 10, 10, 300, 400);
	outputWindow.setVisible(true);
    }

    void saveButton_actionPerformed(ActionEvent event) {
	ExampleFileFilter hfFilter = new ExampleFileFilter("hf","HF HyperFun Data Files");
	ExampleFileFilter wrlFilter = new ExampleFileFilter("wrl","VRML Data Files");
	ExampleFileFilter stlFilter = new ExampleFileFilter("stl","STL Data Files");

	File hf = new File("files/output.hf");
	File wrl = new File("files/output.wrl");
	File stl = new File("files/output.stl");

	if (isHfSaveButtonPushed == true) {
	    fileChooser.addChoosableFileFilter(hfFilter);
	    fileChooser.setSelectedFile(hf);
	    isHfSaveButtonPushed = false;
	    //hfFilter.setExtensionListInDescription(true);
	}

	if (isWrlSaveButtonPushed == true) {
	    fileChooser.addChoosableFileFilter(wrlFilter);
	    fileChooser.setSelectedFile(wrl);
	    isWrlSaveButtonPushed = false;
	}
	if (isStlSaveButtonPushed == true) {
	    fileChooser.addChoosableFileFilter(stlFilter);
	    fileChooser.setSelectedFile(stl);
	    isStlSaveButtonPushed = false;
	}

	fileChooser.setAcceptAllFileFilterUsed(false);
	int returnVal = fileChooser.showSaveDialog(this);

	try {
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		FileWriter fw = new FileWriter(fileChooser.getSelectedFile());
		fw.write(outputtxt);
		fw.close();
	    }
	} catch (Exception ex) {
	    ex.printStackTrace();
	}
    }

    void openButton_actionPerformed(ActionEvent e) {
	ExampleFileFilter hfFilter = new ExampleFileFilter("hf","HF HyperFun Data Files");
	fileChooser.addChoosableFileFilter(hfFilter);
	fileChooser.setAcceptAllFileFilterUsed(false);
	int returnVal = fileChooser.showOpenDialog(this);
	try {
	    if (returnVal == JFileChooser.APPROVE_OPTION) {
		BufferedReader buf = new BufferedReader(new FileReader (fileChooser.getSelectedFile()));
		String tmp, str="";
		while ((tmp = buf.readLine()) != null) {
		    str = str + tmp;
		    str = str + "\n";
		}

		jTextArea1.setText(str);

		jTextArea1.setCaretPosition(0);
		buf.close();
	    }
	} catch (Exception e2) {
	}
    }

    void jpnButton_actionPerformed(ActionEvent event) {
	// Icon on the tab menu
	ImageIcon icon = null;
	TitledBorder title;

	if (toEnglish == true) {
	    jpnButton.setText("Japanese");
	    tabbedPane.addTab("Model Space", icon, panel1, "You can change your model space.");
	    title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Enter the bounding box values.");
	    panel_bdbox.setBorder(title);
	    resetBdBoxButton.setText("reset");
	    tabbedPane.addTab("Grid", icon, panel2, "You can change the grid size.");
	    title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Enter the number of grid lines.");
	    panel_grid.setBorder(title);
	    resetGridButton.setText("reset");
	    tabbedPane.addTab("Color", icon, panel3, "You can change your object color.");
	    title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "Enter the number of grid lines.");
	    panelColor.setBorder(title);
	    selectButton.setText("select");
	    tabbedPane.addTab("Library", icon, panel4, "You can reference HF library.");
	    tabbedPane.addTab("Save", icon, panel5, "You can save the file.");
	    hfSaveButton.setText("Save as HF");
	    hfOpenButton.setText("Open an HF file");
	    wrlSaveButton.setText("Save as VRML");
	    stlSaveButton.setText("Save as STL");
	    povraySaveButton.setText("Save as povray mesh");
	}
	if (toJapanese == true) {
	    jpnButton.setText("  English ");
	    tabbedPane.addTab("\u30e2\u30c7\u30eb\u30b9\u30da\u30fc\u30b9", icon, panel1, "\u30d0\u30a6\u30f3\u30c7\u30a3\u30f3\u30b0\u30dc\u30c3\u30af\u30b9\u306e\u5927\u304d\u3055\u3092\u6c7a\u3081\u307e\u3059\u3002");
	    title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "\u30d0\u30a6\u30f3\u30c7\u30a3\u30f3\u30b0\u30dc\u30c3\u30af\u30b9\u306e\u5024\u3092\u5165\u529b");
	    panel_bdbox.setBorder(title);
	    resetBdBoxButton.setText("\u30ea\u30bb\u30c3\u30c8");
	    tabbedPane.addTab("\u30b0\u30ea\u30c3\u30c9", icon, panel2, "\u30b0\u30ea\u30c3\u30c9\u306e\u5927\u304d\u3055\u3092\u6c7a\u3081\u307e\u3059\u3002");
	    title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "\u30d0\u30a6\u30f3\u30c7\u30a3\u30f3\u30b0\u30dc\u30c3\u30af\u30b9\u306e\u5206\u5272\u6570\u3092\u5165\u529b");
	    panel_grid.setBorder(title);
	    resetGridButton.setText("\u30ea\u30bb\u30c3\u30c8");
	    tabbedPane.addTab(" \u8272 ", icon, panel3, "\u30aa\u30d6\u30b8\u30a7\u30af\u30c8\u306e\u8272\u3092\u5909\u3048\u3089\u308c\u307e\u3059\u3002");
	    title = BorderFactory.createTitledBorder( BorderFactory.createEtchedBorder(EtchedBorder.LOWERED), "\u30aa\u30d6\u30b8\u30a7\u30af\u30c8\u306e\u8272\u3092\u9078\u629e");
	    panelColor.setBorder(title);
	    selectButton.setText("\u9078\u629e");
	    tabbedPane.addTab("\u30e9\u30a4\u30d6\u30e9\u30ea", icon, panel4, "\u30e9\u30a4\u30d6\u30e9\u30ea\u306e\u53c2\u7167\u3002");
	    tabbedPane.addTab("\u4fdd\u5b58", icon, panel5, "\u4fdd\u5b58\u3067\u304d\u307e\u3059\u3002");
	    hfSaveButton.setText("HF\u3067\u4fdd\u5b58");
	    wrlSaveButton.setText("VRML\u3067\u4fdd\u5b58");
	    stlSaveButton.setText("STL\u3067\u4fdd\u5b58");
	}

	if (toEnglish == true) {
	    toEnglish = false;
	    isEnglish = true;
	    isJapanese = false;
	}
	if (toJapanese == true) {
	    toJapanese = false;
	    isEnglish = false;
	    isJapanese = true;
	}
    }
}
