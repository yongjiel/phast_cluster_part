import java.awt.*;
import java.awt.image.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import javax.imageio.*;

public class Plotter {
	public static void main(String args[]) {
		if (args.length > 0) {
			Plotter plot = new Plotter();
			String line;
			String[] parts;
			int id = 0, level = 0;
			try {
				FileReader input = new FileReader(args[0]);
				BufferedReader bufRead = new BufferedReader(input);
				while ((line = bufRead.readLine()) != null){
					if (line.startsWith(">")){
						parts = line.split(",");
						id = plot.addPlasmid(new Integer(parts[parts.length-1].trim()), level++, parts[0].substring(1));
					}
					if (line.startsWith("section")){
						parts = line.split("\\s+");
						plot.addGene(id, new Integer(parts[1]), new Integer(parts[2]), parts[3], parts[4]);
					}
					if (line.startsWith("region")){
						parts = line.split("\\s+");
						plot.addRegion(id, new Integer(parts[2]), new Integer(parts[3]), parts[4]);
					}
				}
			} catch (Exception e) {// Catch exception if any
				System.err.println("Error: " + e.getMessage());
			}
			if (args.length >= 3 && args[2].equals("-a"))
				plot.plasmids.changeReservedSymbol();
			plot.resize();
			plot.plot();
			if (args.length > 1)
				plot.save(args[1]);
			else
				plot.save("figure.png");
		}
	}
	
	public static final int[] SIZE = {400, 100, 50};
	private BufferedImage plot;
	private Graphics2D graph;
	private int width;
	private int height;
	private Plasmid plasmids;
	private int numPlasmid = 0;
	public Plotter() {
		this.width = 1500;
		this.height = 1000;
		this.plot = new BufferedImage(this.width, this.height,
				BufferedImage.TYPE_INT_RGB);
		this.graph = this.plot.createGraphics();
		this.graph.setColor(Color.WHITE);
		this.graph.fillRect(0, 0, this.width, this.height);
	}
	/** resize the figure to include detailed regions */
	public void resize(){
		this.width = 1700;
		this.height = 1000;
		int numr = 0;
		Plasmid p = this.plasmids;
		while (p != null){
			numr += p.numOfRegions();
			p = p.next;
		}
		this.height += 150 * numr;
		this.plot = new BufferedImage(this.width, this.height,
				BufferedImage.TYPE_INT_RGB);
		this.graph = this.plot.createGraphics();
		this.graph.setColor(Color.WHITE);
		this.graph.fillRect(0, 0, this.width, this.height);
	}
	public Graphics2D getGraphics2D() {
		return this.graph;
	}

	public int getWidth() {
		return this.width;
	}

	public int getHeight() {
		return this.height;
	}

	public void save(String filename) {
		try {
			File outputfile = new File(filename);
			ImageIO.write(this.plot, "png", outputfile);
		} catch (IOException e) {
			System.err.println("cannot write to file: " + filename);
		}

	}
	
	public int addPlasmid(int bp, int level, String name){
		this.numPlasmid++;
		Plasmid np = new Plasmid(this, this.numPlasmid, bp, SIZE[level], 1000, 500);
		np.setName(name);
		if (this.plasmids == null)
			this.plasmids = np;
		else{
			Plasmid last = this.plasmids;
			while (last.next() != null)
				last = last.next();
			last.setNext(np);
		}
		return this.numPlasmid;
	}
	
	public int addGene(int plasmidId, int from, int to, String code, String name){
		Plasmid p = this.plasmids;
		int frame;
		boolean direct;
		if (code.startsWith("+"))
			direct = Gene.OUT;
		else if (code.startsWith("-"))
			direct = Gene.IN;
		else
			return -1;
		frame = new Integer(code.substring(1));
		for (;p != null;p=p.next())
			if (p.getId() == plasmidId)
				p.addGene(from, to, frame, direct, name);
		return 0;
	}
	public void addRegion(int plasmidId, int from, int to, String type){
		Plasmid p = this.plasmids;
		boolean t;
		if (type.equals("true"))
			t = true;
		else
			t = false;
		for (;p != null;p=p.next())
			if (p.getId() == plasmidId)
				p.addRegion(from, to, t);
	}
	public void plot(){
		for(Plasmid p = this.plasmids; p!=null; p=p.next()){
			p.show();
		}
	}
	public void tplot() {
		
	}
}
