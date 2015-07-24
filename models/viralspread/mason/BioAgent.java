package sim.app.virus3;

import sim.engine.*;
import sim.util.Double2D;
import sim.portrayal.*;
import java.awt.geom.*;

// Based on Agent
public abstract class BioAgent extends sim.portrayal.SimplePortrayal2D implements Steppable
    {
    private static final long serialVersionUID = 1;
    public String id;
    public Double2D agentLocation; 
    public int intID = -1;
    
    public BioAgent( SimpleVirusDemo agent, String id, Double2D location )
        {
        this.id = id;
        this.agentLocation = location;
        }

    double distanceSquared( final Double2D loc1, Double2D loc2 )
        {
        return( (loc1.x-loc2.x)*(loc1.x-loc2.x)+(loc1.y-loc2.y)*(loc1.y-loc2.y) );
        }

    // Returns "Cell" or "Virion"
    public abstract String getType();  

    public boolean hitObject(Object object, DrawInfo2D info)
        {
        double diamx = info.draw.width*SimpleVirusDemo.CELL_DIAMETER;
        double diamy = info.draw.height*SimpleVirusDemo.CELL_DIAMETER;

        Ellipse2D.Double ellipse = new Ellipse2D.Double( (int)(info.draw.x-diamx/2),(int)(info.draw.y-diamy/2),(int)(diamx),(int)(diamy) );
        return ( ellipse.intersects( info.clip.x, info.clip.y, info.clip.width, info.clip.height ) );
        }
    }