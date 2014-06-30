
public class Region {
	public int from;
	public int to;
	public int num;
	public Region next;
	public boolean type;
	public Region(int from, int to, int num, boolean type){
		this.from = from;
		this.to = to;
		this.num = num;
		this.type = type;
	}
}
