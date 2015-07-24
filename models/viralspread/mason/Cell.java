package sim.app.virus3;

import sim.util.*;
import sim.engine.*;
import java.awt.*;
import sim.portrayal.*;

// Based on Human
public class Cell extends BioAgent
    {
    private static final long serialVersionUID = 1;
    
    // If the state is true the cell is infected
    protected boolean infected = false;
    // This method checks the status of the cell
    public final boolean isInfected() { return infected; }
    // This method may be used to change the status of the cell
    public final void setInfected( boolean b ) { infected = b; }
    
    // Main class
    public Cell( SimpleVirusDemo agent, String id, Double2D location ) 
        {
        super( agent, id, location );
        try
            {
            // This returns an integer starting at the fifth location in the string id
            // Changed 5 to 4 (length of word Human vs Cell)
            // The id is set in SimpleVirusDemo when a new agent is made
            intID = Integer.parseInt( id.substring(4) ); // "Cell"
            }
        catch( Exception e )
            {
            throw new RuntimeException(e);
            }
        }
    
    public void step( final SimState state )
        {
        }
    
    // R G B
    //protected Color healthyColor = new Color(192,128,128);
    protected Color healthyColor = new Color(0,255,0);
    protected Color infectedColor = new Color(0,0,255);
    public final void draw(Object object, Graphics2D graphics, DrawInfo2D info)
        {
        double diamx = info.draw.width*SimpleVirusDemo.CELL_DIAMETER;
        double diamy = info.draw.height*SimpleVirusDemo.CELL_DIAMETER;
    
        if (isInfected())
            graphics.setColor( infectedColor );
        else
            graphics.setColor ( healthyColor );
        
        graphics.fillOval((int)(info.draw.x-diamx/2),(int)(info.draw.y-diamy/2),(int)(diamx),(int)(diamy));
        }

    public String getType()
        {
        if( isInfected() )
            return "Infected Cell";
        else
            return "Healthy Cell";
        }
    }
