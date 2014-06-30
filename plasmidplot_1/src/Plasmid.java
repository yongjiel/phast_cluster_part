import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.geom.*;
import java.util.Hashtable;
import java.util.Random;

public class Plasmid {
	protected static final Color BACTERIA = Color.black;
	protected static final Color HYPOTHETICAL = Color.darkGray;
	protected Random numGen = new Random(83);
	protected Plotter parent;
	protected Plasmid next;
	protected final int id;
	protected String name;
	protected int size;
	protected float rad;
	protected int centx;
	protected int centy;
	protected int majorscale = 12;
	protected int minorscale = 60;
	protected Font labelfont = new Font("Arial", Font.BOLD, 12);
	protected static final BasicStroke framestrock = new BasicStroke(1f);
	protected static final BasicStroke majorstrock = new BasicStroke(3f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	protected static final BasicStroke minorstrock = new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
	protected Gene genes;
	protected int numgene = 0;
	protected Region regions;
	protected int regionspace = 130;
	protected int numColor = 0;
	protected Hashtable<Color, Integer> uniquecolors = new Hashtable<Color, Integer>(128);
	protected Hashtable<String, Color> uniquegenes = new Hashtable<String, Color>(128);
	protected String reservedSymbol1 = "b";
	protected String reservedSymbol2 = "h";
	public Plasmid(Plotter parent, int id, int size, float rad, int x, int y){
		this.centx = x;
		this.centy = y;
		this.parent = parent;
		this.id = id;
		this.size = size;
		this.rad = rad;
		this.uniquecolors.put(BACTERIA, 0);
//		this.uniquegenes.put("Bacterial_gene", BACTERIA);
		this.uniquecolors.put(HYPOTHETICAL, 1);
//		this.uniquegenes.put("Hypothetical_gene", HYPOTHETICAL);
	}
	private String getNumLabel(float number){
		if (number >= 1000000f)
			return String.format("%.1fm", number/1000000f);
		else if (number >= 1000f){
			return String.format("%.1fk", number/1000f);
//			return new Integer((int)number/1000).toString()+"k";
		}
		return new Integer((int)number).toString();
	}
	public void setName(String name){
		this.name = name;
	}
	public int numOfRegions(){
		int num = 0;
		Region r = this.regions;
		while (r != null){
			num++;
			r = r.next;
		}
		return num;
	}
	
	/**
	 * draw figure
	 * 
	 */
	public void show(){
		Graphics2D graph = this.parent.getGraphics2D();
		float level, angle, from, to;
		GeneralPath gp;
		Color color;
		// draw frame circle
		graph.setColor(Color.black);
		graph.setStroke(framestrock);
		this.centy = this.getRingDepth();
		graph.drawArc((int)(this.centx - this.rad),
				(int)(this.centy - this.rad),
				(int)(2*this.rad), (int)(2*this.rad),
				0, 360);
		graph.drawString(this.size+" bp", this.centx-30, this.centy+5);
		graph.drawString(this.name, this.centx-180, this.centy+20);
		// draw major scales
		graph.setStroke(majorstrock);
		graph.setFont(this.labelfont);
		for (float i = 0; i < majorscale; i ++){
			graph.drawLine((int)(this.getX(rad+40, 2*(float)Math.PI*i/majorscale)),
					(int)(this.getY(rad+40, 2*(float)Math.PI*i/majorscale)),
					(int)(this.getX(rad+40, 2*(float)Math.PI*i/majorscale)),
							(int)(this.getY(rad+40, 2*(float)Math.PI*i/majorscale)));
			graph.drawString(getNumLabel(this.size * i/majorscale),
					(int)(this.getX(rad+55, -2*(float)Math.PI*i/majorscale+(float)Math.PI/2)) - 10, 
					(int)(this.getY(rad+55, -2*(float)Math.PI*i/majorscale+(float)Math.PI/2)) + 10);
		}
		graph.setStroke(minorstrock);
		for (float i = 0; i < minorscale; i ++){
			graph.drawLine((int)((rad+40)*Math.sin(2*Math.PI*i/minorscale) + this.centx), 
					(int)((rad+40)*Math.cos(2*Math.PI*i/minorscale) + this.centy),
					(int)((rad+40)*Math.sin(2*Math.PI*i/minorscale) + this.centx), 
					(int)((rad+40)*Math.cos(2*Math.PI*i/minorscale)) + this.centy);
		}
		// mark gene
		int yscale =70;
		for (Gene g = this.genes; g != null; g = g.next()){
			gp = new GeneralPath();
			graph.setStroke(Plasmid.minorstrock);
			if (g.getDiection()){ // OUT SIZE
				if (g.getFrame() == 0)
					color = Plasmid.BACTERIA;
				else if (g.getFrame() == 1)
					color = Plasmid.HYPOTHETICAL;
				else if (this.uniquegenes.get(g.getName()) != null)
					color = this.uniquegenes.get(g.getName());
				else 
					color = this.getNewColor();
				graph.setPaint(color);
				level = this.rad + 8;
				from = this.getAngle(g.from());
				to = this.getAngle(g.to());
				if (from < 0)
					from += Math.PI*2;
				if (to < 0)
					to += Math.PI*2;
				gp.moveTo(this.getX(level, from), this.getY(level, from));
				for (angle = from; angle > to; angle -= 1f/360f)
					gp.lineTo(this.getX(level, angle), this.getY(level, angle));
				gp.lineTo(this.getX(level, to), this.getY(level, to));
				level += 20;
				gp.lineTo(this.getX(level, to), this.getY(level, to));
				for (angle = to; angle < from; angle += 1f/360f)
					gp.lineTo(this.getX(level, angle), this.getY(level, angle));
				gp.lineTo(this.getX(level, from), this.getY(level, from));
				level -= 20;
				gp.lineTo(this.getX(level, from), this.getY(level, from));
				graph.fill(gp);
			}
			else{
				if (g.getFrame() == 0)
					color = Plasmid.BACTERIA;
				else if (g.getFrame() == 1)
					color = Plasmid.HYPOTHETICAL;
				else if (this.uniquegenes.get(g.getName()) != null)
					color = this.uniquegenes.get(g.getName());
				else 
					color = this.getNewColor();
				graph.setPaint(color);
				level = this.rad - 8;
				from = this.getAngle(g.from());
				to = this.getAngle(g.to());
				if (from < 0)
					from += Math.PI*2;
				if (to < 0)
					to += Math.PI*2;
				gp.moveTo(this.getX(level, from), this.getY(level, from));
				for (angle = from; angle > to; angle -= 1f/360f)
					gp.lineTo(this.getX(level, angle), this.getY(level, angle));
				gp.lineTo(this.getX(level, to), this.getY(level, to));
				level -= 20;
				gp.lineTo(this.getX(level, to), this.getY(level, to));
				for (angle = to; angle < from; angle += 1f/360f)
					gp.lineTo(this.getX(level, angle), this.getY(level, angle));
				gp.lineTo(this.getX(level, from), this.getY(level, from));
				level += 20;
				gp.lineTo(this.getX(level, from), this.getY(level, from));
				graph.fill(gp);			
			}
			// legend
			if (uniquegenes.get(g.getName()) != null)
				continue;
			uniquegenes.put(g.getName(), color);
			graph.setStroke(Plasmid.majorstrock);
			if (g.getFrame() == 0){
				graph.drawLine(30, 40, 50, 40);
				graph.setStroke(Plasmid.framestrock);
				graph.drawString(this.reservedSymbol1, 60, 45);
				graph.drawString(g.getName(), 90, 45);
			}
			else if (g.getFrame() == 1){
				graph.drawLine(30, 55, 50, 55);
				graph.setStroke(Plasmid.framestrock);
				graph.drawString(this.reservedSymbol2, 60, 60);
				graph.drawString(g.getName(), 90, 60);
			}
			else{
				graph.drawLine(30, yscale, 50, yscale);
				graph.setStroke(Plasmid.framestrock);
				graph.drawString(this.uniquecolors.get(color).toString(), 60, yscale+5);
				graph.drawString(g.getName(), 90, yscale+5);
				yscale += 15;
			}
		}
		// regions
		Region r;
		AffineTransform old;
		AffineTransform current;
		int depth = 80;
		int left = this.centx - (int)this.rad - 50;
		int right = this.centx + (int)this.rad + 50;
//		graph.drawString("Note: regions may have different scales. " +
//				"Band size should be estimated by base pair per unit legend to the right.", left-50, depth-110);
//		graph.drawString("Note: this figure shows both true and defective regions. " +
//				"You can hide defective regions by a parameter when run the plotter.", left-50, depth-90);
		for (r = this.regions; r != null; r = r.next){
			graph.setStroke(Plasmid.minorstrock);
			graph.setPaint(Color.black);
			//draw arc
			angle = this.getAngle(r.from);
			int lb = (int)(this.getAngle(r.to)*180/Math.PI);
			int ub = (int)(this.getAngle(r.from)*180/Math.PI);
			graph.drawArc((int)(this.centx - this.rad+80),
				(int)(this.centy - this.rad+80),
				(int)(2*(this.rad-80)), (int)(2*(this.rad-80)),
				lb, Math.max(1, ub - lb));
			graph.drawLine((int)this.getX(this.rad-85, angle), (int)this.getY(this.rad-85, angle),
					(int)this.getX(this.rad-80, angle), (int)this.getY(this.rad-80, angle));
			old = graph.getTransform();
			current = AffineTransform.getRotateInstance(Math.PI-angle, this.getX(this.rad-90, angle), this.getY(this.rad-90, angle));
			graph.transform(current);
			graph.drawString("Region "+new Integer(r.num).toString(), this.getX(this.rad-90, angle), this.getY(this.rad-90, angle));
			graph.setTransform(old);
			// draw line
			graph.drawLine(left, depth, right, depth);
			// draw scale
			for (int i = left; i <= right; i+= (right-left)/10)
				graph.drawLine(i, depth+5, i, depth-5);
			for (int i = left; i <= right; i+= (right-left)/20)
				graph.drawLine(i, depth+2, i, depth-2);
			/*
			graph.drawLine(right+100, depth-30, right+100+(right-left)/20, depth-30);
			graph.drawLine(right+100, depth-30, right+100, depth-32);
			graph.drawLine(right+100+(right-left)/20, depth-30, right+100+(right-left)/20, depth-32);
			graph.drawString((r.to - r.from + 1)/20+" bp", right+100, depth-12);
			graph.drawLine(right+100, depth, right+100+(right-left)/10, depth);
			graph.drawLine(right+100, depth, right+100, depth-5);
			graph.drawLine(right+100+(right-left)/20, depth, right+100+(right-left)/20, depth-2);
			graph.drawLine(right+100+(right-left)/10, depth, right+100+(right-left)/10, depth-5);
			graph.drawString((r.to - r.from + 1)/10+" bp", right+100, depth+18);
			*/
			graph.drawLine(right+100, depth, right+100+(right-left)/10, depth);
			graph.drawLine(right+100, depth, right+100, depth-5);
			graph.drawLine(right+100+(right-left)/20, depth, right+100+(right-left)/20, depth-2);
			graph.drawLine(right+100+(right-left)/10, depth, right+100+(right-left)/10, depth-5);
			graph.drawString("0 bp", right+100, depth+18);
			graph.drawString((r.to - r.from + 1)/20+" bp", right+100+(right-left)/20, depth-10);
			graph.drawString((r.to - r.from + 1)/10+" bp", right+100+(right-left)/10, depth+18);
			// draw number	
			if (r.type){
				graph.setPaint(Color.green.darker());
				graph.drawString("true", left-130, depth+5);
			}
			else{
				graph.setPaint(Color.red);
				graph.drawString("defective", left-130, depth+5);
			}
			graph.drawString("Region "+r.num, left-190, depth+5);
			graph.setPaint(Color.black);
			graph.drawString(new Integer(r.from).toString(), left-60, depth+5);
			graph.drawString(new Integer(r.to).toString(), right+10, depth+5);
			// draw arrows
			/*
			gp = new GeneralPath();
			gp.moveTo(left-110, depth-25);
			gp.lineTo(left-40, depth-25);
			gp.lineTo(left-40, depth-30);
			gp.lineTo(left-20, depth-20);
			gp.lineTo(left-40, depth-10);
			gp.lineTo(left-40, depth-15);
			gp.lineTo(left-110, depth-15);
			gp.lineTo(left-110, depth-25);
			graph.draw(gp);
			gp = new GeneralPath();
			gp.moveTo(left-20, depth+25);
			gp.lineTo(left-90, depth+25);
			gp.lineTo(left-90, depth+30);
			gp.lineTo(left-110, depth+20);
			gp.lineTo(left-90, depth+10);
			gp.lineTo(left-90, depth+15);
			gp.lineTo(left-20, depth+15);
			gp.lineTo(left-20, depth+25);
			graph.draw(gp);
			*/
			// draw band
			String label;
			int bandnum = 0;
			Rectangle2D[] bands = new Rectangle2D[512];
			Color[] bandcolor = new Color[512];
			float h = 0f;
			for (Gene g = this.genes; g != null; g = g.next()){
				if (g.from() >= r.from && g.to() <= r.to){
					gp = new GeneralPath();
					color = this.uniquegenes.get(g.getName());
					if (g.getFrame() == 0)
						label = "b";
					else if (g.getFrame() == 1)
						label = "h";
					else
						label = this.uniquecolors.get(color).toString();
					graph.setPaint(color);
					bandcolor[bandnum] = new Color(color.getRed(), color.getGreen(), color.getBlue());
					if (bandnum > 0 &&
							bands[bandnum-1].getX()+bands[bandnum-1].getWidth() > this.getLineX(g.from(), r, left, right)
							&& ((bands[bandnum-1].getY() < depth && g.getDiection()) ||
									(bands[bandnum-1].getY() > depth && !g.getDiection()))){
						if (h == 0f)
							h = 5f;
						else
							h = 0f;
					}
					else
						h = 0f;
					if (g.getDiection()){
						level = depth - 8;
						bands[bandnum] = new Rectangle2D.Float(this.getLineX(g.from(), r, left, right), level - 20 + h,
								this.getLineX(g.to(), r, left, right)- this.getLineX(g.from(), r, left, right),
								20-h);
//						graph.fillRect((int)this.getLineX(g.from(), r, left, right), (int)level - 20,
//								(int)this.getLineX(g.to(), r, left, right)- (int)this.getLineX(g.from(), r, left, right),
//								20);
						graph.drawString(label, (int)this.getLineX(g.from(), r, left, right), (int)level - 25);
						graph.setPaint(Color.black);
						graph.drawLine((int)this.getLineX(g.from(), r, left, right), (int)(level - 20 + h),
								(int)this.getLineX(g.from(), r, left, right), depth);
						graph.drawLine((int)this.getLineX(g.to(), r, left, right), (int)(level - 20 + h),
								(int)this.getLineX(g.to(), r, left, right), depth);
					}
					else{
						level = depth + 8;
						bands[bandnum] = new Rectangle2D.Float(this.getLineX(g.from(), r, left, right), level,
								this.getLineX(g.to(), r, left, right)- this.getLineX(g.from(), r, left, right),
								20-h);
//						graph.fillRect((int)this.getLineX(g.from(), r, left, right), (int)level,
//								(int)this.getLineX(g.to(), r, left, right)- (int)this.getLineX(g.from(), r, left, right),
//								20);
						graph.drawString(label, (int)this.getLineX(g.from(), r, left, right), (int)level + 35);
						graph.setPaint(Color.black);
						graph.drawLine((int)this.getLineX(g.from(), r, left, right), (int)(level + 20 - h),
								(int)this.getLineX(g.from(), r, left, right), depth);
						graph.drawLine((int)this.getLineX(g.to(), r, left, right), (int)(level + 20 - h),
								(int)this.getLineX(g.to(), r, left, right), depth);
					}
					bandnum++;
				}
			}
			// paint bands based on their sizes (paint big bands first)
			// sort band arrays
			for (int x = 0; x < bandnum; x++)
				for (int y = x+1; y < bandnum; y++)
					if (bands[x].getHeight() < bands[y].getHeight()){
						Rectangle2D b = bands[x];
						bands[x] = bands[y];
						bands[y] = b;
						Color c = bandcolor[x];
						bandcolor[x] = bandcolor[y];
						bandcolor[y] = c;
					}
			// paint bands
			for (int i = 0; i < bandnum; i++){
				graph.setPaint(bandcolor[i]);
				graph.fillRect((int)bands[i].getX(), (int)bands[i].getY(), (int)bands[i].getWidth(), (int)bands[i].getHeight());
	//			graph.fill(bands[i]);
				graph.setPaint(Color.black);
	//			graph.draw(bands[i]);
				graph.drawRect((int)bands[i].getX(), (int)bands[i].getY(), (int)bands[i].getWidth(), (int)bands[i].getHeight());
			}
			depth += this.regionspace;
		}
		
	}
	public void addGene(int from, int to, int frame, boolean direct, String name){
		if (this.numgene == 0)
			this.genes = new Gene(from, to, frame, direct, name);
		else{
			Gene last = this.genes;
			while (last.next() != null)
				last = last.next();
			last.setNext(new Gene(from, to, frame, direct, name));
		}
		this.numgene++;
	}
	public Plasmid next(){
		return this.next;
	}
	public void setNext(Plasmid next){
		this.next = next;
	}
	public int getId(){
		return this.id;
	}
	protected float getAngle(int index){
		return (float)(Math.PI/2 - 2 * Math.PI * index / this.size);
	}
	protected float getX(float length, float angle){
		return this.centx + (float)(length * Math.cos(angle));
	}
	protected float getY(float length, float angle){
		return this.centy - (float)(length * Math.sin(angle));
	}
	protected float getLineX(int value, Region r, int left, int right){
		return (value - r.from)/(float)(r.to - r.from)*(right - left)+left;
	}
	protected Color getNewColor(){
		Color nc = new Color(numGen.nextInt(256), numGen.nextInt(256), numGen.nextInt(256));
		while (this.uniquecolors.get(nc) != null || nc.equals(Color.white))
			nc = new Color(numGen.nextInt(256), numGen.nextInt(256), numGen.nextInt(256));
		this.uniquecolors.put(nc, ++this.numColor);
		return nc;
	}
	public void addRegion(int from, int to, boolean type){
		Region r = this.regions;
		int num = 1;
		if (r == null)
			this.regions = new Region(from, to, num, type);
		else{
			num = 2;
			while (r.next != null){
				r = r.next;
				num++;
			}
			r.next = new Region(from, to, num, type);
		}
	}
	protected int getRingDepth(){
		int start = this.centy;
		Region r = this.regions;
		while (r != null){
			start += this.regionspace;
			r = r.next;
		}
		return start+100;
	}
	/**
	 * Jule 16 2010, a flag can change legend symbol 'h', 'b' to 'u', 'a'
	 */
	protected void changeReservedSymbol(){
		this.reservedSymbol1 = "u";
		this.reservedSymbol2 = "a";
	}
}
