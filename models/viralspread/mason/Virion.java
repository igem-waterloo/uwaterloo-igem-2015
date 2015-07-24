package sim.app.virus3;

import sim.util.*;
import sim.engine.*;
import java.awt.*;
import java.util.Random;
import sim.portrayal.*;


public class Virion extends BioAgent
    {
    private static final long serialVersionUID = 1;
    
    protected boolean inNucleus = false;
    public final boolean getInNucleus() { return inNucleus; }
    public final void setInNucleus( final boolean b ) { inNucleus = b; }
    
    protected boolean greedy = true;
    public final boolean getIsGreedy() { return greedy; }
    public final void setIsGreedy( final boolean b ) { greedy = b; }

    public Virion( SimpleVirusDemo agent, String id, Double2D location ) 
        {
        super( agent, id, location );
        try
            {
            intID = Integer.parseInt( id.substring(6) ); // "Cell"
            }
        catch( Exception e )
            {
            throw new RuntimeException(e);
            }
        agent.environment.setObjectLocation(this, location);
        //agent.virionCount++;
        }
    
    Double2D desiredLocation = null;
    Double2D suggestedLocation = null;
    Double2D finalLocation = null;
    int steps = 0;

    public void step( final SimState state )
        {
        SimpleVirusDemo hb = (SimpleVirusDemo)state;

        desiredLocation = null;
        double distance2DesiredLocation = 1e30;
        
        Bag mysteriousObjects = hb.environment.getNeighborsWithinDistance(agentLocation, 2.0 * SimpleVirusDemo.INFECTION_DISTANCE);
        if( mysteriousObjects != null )
            {
            for( int i = 0 ; i < mysteriousObjects.numObjs ; i++ )
                {
                if( mysteriousObjects.objs[i] != null &&
                    mysteriousObjects.objs[i] != this )
                    {
                    // if agent is not cell, wasted time....
                    if( ! (((BioAgent)mysteriousObjects.objs[i]) instanceof Cell ))
                        continue;
                    Cell ta = (Cell)(mysteriousObjects.objs[i]);
                    // if agent is already infected, wasted time....
                    if( ta.isInfected() )
                        continue;
                    else if( hb.withinInfectionDistance( this, agentLocation, ta, ta.agentLocation ) )
                        {
                        ta.setInfected( true );
                        this.setInNucleus( true );
                        //hb.addNew( true );
                        finalLocation = ta.agentLocation;
                        Double2D loc = null;
                        loc = ta.agentLocation;
                        addVirions( hb, loc);
                        //addVirion( hb, loc);
                        //agent.schedule.scheduleRepeating(agent);
                        }
                    else
                        {
                            if( getIsGreedy() )
                            {
                            double tmpDist = distanceSquared( agentLocation, ta.agentLocation );
                            if( tmpDist <  distance2DesiredLocation )
                                {
                                desiredLocation = ta.agentLocation;
                                distance2DesiredLocation = tmpDist;
                                }
                            }
                        }
                    }
                }
            }
        steps--;
        if( desiredLocation == null || !getIsGreedy() )
            {
            if(  steps <= 0 )
                {
                suggestedLocation = new Double2D((state.random.nextDouble()-0.5)*((SimpleVirusDemo.XMAX-SimpleVirusDemo.XMIN)/5-SimpleVirusDemo.VIR_DIAMETER) +
                    //VirusInfectionDemo.XMIN
                    agentLocation.x 
                    //+VirusInfectionDemo.DIAMETER/2
                    ,
                    (state.random.nextDouble()-0.5)*((SimpleVirusDemo.YMAX-SimpleVirusDemo.YMIN)/5-SimpleVirusDemo.VIR_DIAMETER) +
                    agentLocation.y
                    //VirusInfectionDemo.YMIN
                    //+VirusInfectionDemo.DIAMETER/2
                    );
                steps = 100;
                }
            desiredLocation = suggestedLocation;
            }
        double dx = desiredLocation.x - agentLocation.x;
        double dy = desiredLocation.y - agentLocation.y;
        
                {
                double temp = 0.5 * /*Strict*/Math.sqrt(dx*dx+dy*dy);
                if( temp < 1 )
                    {
                    steps = 0;
                    }
                else
                    {
                    dx /= temp;
                    dy /= temp;
                    }
                }

        if ( getInNucleus() )
            {
            //agentLocation = new Double2D(agentLocation.x, agentLocation.y);
            hb.environment.setObjectLocation(this,finalLocation);
            }
        else
            {
            agentLocation = new Double2D(agentLocation.x + dx, agentLocation.y + dy);
            hb.environment.setObjectLocation(this,agentLocation);
            }
        }

    protected Color virionColor = new Color(255,0,0);
    protected Color virionNucleusColor = new Color(0,0,0);
    public final void draw(Object object, Graphics2D graphics, DrawInfo2D info)
        {
        double diamx = info.draw.width*SimpleVirusDemo.VIR_DIAMETER;
        double diamy = info.draw.height*SimpleVirusDemo.VIR_DIAMETER;
        
        if ( getInNucleus() )
            graphics.setColor( virionNucleusColor );
        else
            {
            graphics.setColor( virionColor );
            graphics.fillOval((int)(info.draw.x-diamx/2),(int)(info.draw.y-diamy/2),(int)(diamx),(int)(diamy));
            }
        }
    
    public void addVirions( final SimpleVirusDemo hb, Double2D loc )
    {
        BioAgent agent1 = null;
        agent1 = new Virion( hb, "Virion" + hb.virionCount, loc);
        hb.virionCount++;
        hb.schedule.scheduleRepeating(agent1);
        
        double randNum = Math.random();
        if ( randNum < 0.4 )
            {
            BioAgent agent2 = null;
            agent2 = new Virion( hb, "Virion" + hb.virionCount, loc);
            hb.virionCount++;
            hb.schedule.scheduleRepeating(agent2);
            }
    }
    
    public String getType() { return "Virion"; }

    }
