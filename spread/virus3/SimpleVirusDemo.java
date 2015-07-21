package sim.app.virus3;

import sim.field.continuous.*;
import sim.engine.*;
import sim.util.*;

// Based on SimpleVirusDemo
public class SimpleVirusDemo extends SimState
    {
    private static final long serialVersionUID = 1;

    public static final double XMIN = 0;
    public static final double XMAX = 400;
    public static final double YMIN = 0;
    public static final double YMAX = 400;
    public static final double XLEN = XMAX-XMIN;
    public static final double YLEN = YMAX-YMIN;
    
    public static final double CELL_DIAMETER = 10;
    public static final double VIR_DIAMETER = 3;
    
    public static final double INFECTION_DISTANCE = 3;//20;
    public static final double INFECTION_DISTANCE_SQUARED = INFECTION_DISTANCE * INFECTION_DISTANCE;
    
    public static final int NUM_CELLS_X = 40;
    public static final int NUM_CELLS_Y = 40;
    public static final int NUM_CELLS = NUM_CELLS_X * NUM_CELLS_X;
    public static final int NUM_VIRIONS = 1;
    public static final int NUM_INIT = NUM_CELLS + NUM_VIRIONS;
    
    public int virionCount = 0;
    
    // Should a cell be added or not
    protected boolean newAgent = true;
    public final void addNew( boolean b ) { newAgent = b; }
    
    public Continuous2D environment = null;

    /** Creates a VirusInfectionDemo simulation with the given random number seed. */
    public SimpleVirusDemo(long seed)
        {
        super(seed);
        }
    
    // Method to determine if a cell and a virion are close enough for the cell to be infected
    public boolean withinInfectionDistance( final BioAgent agent1, final Double2D a, final BioAgent agent2, final Double2D b )
        {
        return ( (a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y) <= INFECTION_DISTANCE_SQUARED );
        }
    
    // *** Simulation starts here *** //
    
    // start is called immediately before the schedule is iterated
    // We override this method to setup our simulation
    public void start()
        {
        super.start();  // clear out the schedule 

        environment = new Continuous2D(25.0, XLEN, YLEN );
        
        // Add in new virions
        
        // Schedule our agents
        for(int x=0;x<(NUM_INIT+virionCount);x++)
            {
            Double2D loc = null;
            BioAgent agent = null;
            
                if( x < NUM_CELLS )
                    {loc = new Double2D( (x%NUM_CELLS_X)*XLEN/NUM_CELLS_X + XLEN/(2*NUM_CELLS_X),
                                         (x/NUM_CELLS_Y)*YLEN/NUM_CELLS_Y + YLEN/(2*NUM_CELLS_Y));
                    agent = new Cell( this, "Cell"+x, loc );}
                else
                    {loc = new Double2D( random.nextDouble()*XLEN, random.nextDouble()*YLEN);
                    //agent = new Virion( this, "Virion"+(x-NUM_CELLS), loc );}
                    agent = new Virion( this, "Virion"+(x-NUM_CELLS), loc );}
                
            environment.setObjectLocation(agent,loc);
            schedule.scheduleRepeating(agent);
            }
        }
    
    // Used in all top-level simulations
    public static void main(String[] args)
        {
        doLoop(SimpleVirusDemo.class, args);
        System.exit(0);
        }    
    }
