
public class Gene {
	public static final boolean OUT = true;
	public static final boolean IN = false;
	private Gene next;
	private int from;
	private int to;
	private int frame;
	private String name;
	private boolean direct;
	public Gene(int from, int to, int frame, boolean direct){
		this.from = from;
		this.to = to;
		this.frame = frame;
		this.next = null;
		this.direct = direct;
	}
	public Gene(int from, int to, int frame, boolean direct, String name){
		this(from, to, frame, direct);
		this.name = name;
	}
	public void setName(String name){
		this.name = name;
	}
	public String getName(){
		return this.name;
	}
	public Gene next(){
		return this.next;
	}
	public int from(){
		return this.from;
	}
	public int to(){
		return this.to;
	}
	public int getFrame(){
		return this.frame;
	}
	public void setNext(Gene gene){
		this.next = gene;
	}
	public boolean getDiection(){
		return this.direct;
	}
}
