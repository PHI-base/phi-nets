package uk.ac.rothamsted.phinets.util;

import java.util.Arrays;

public class ArrayKey<T> {
	private T[] array = null;
	private int hash;

	public ArrayKey() {
	}
	
	public ArrayKey(boolean sort, T ... array){
		this.array = array;
		sort = true;
		if(sort){
			Arrays.sort(array);	
		}
		this.hash = Arrays.hashCode(this.array);
	}

	public ArrayKey(T ... array){
		this.array = array;
		this.hash = Arrays.hashCode(this.array);
	}
	
	public ArrayKey<T> extend(T ... values){
		T[] arr = Arrays.copyOf(array, array.length+values.length);
		for(int i = 0; i < values.length; i++){
			arr[array.length+i] = values[i];
		}
		ArrayKey<T> key = new ArrayKey<T>(arr);
		return key;
	}

	public T[] getArray() {
		return this.array;
	}

	public void setArray(T[] array) {
		this.array = array;
		this.hash = Arrays.hashCode(this.array);
	}

	public int hashCode() {
		return hash;
	}

	public int size() {
		return this.array.length;
	}

	public boolean equals(Object key) {
		if(this == key){
			return true;
		}
		if ((key instanceof ArrayKey)) {
			return Arrays.equals(this.array, ((ArrayKey) key).array);
		}
		return false;
	}
	
	public void sort(){
		Arrays.sort(array);
		this.hash = Arrays.hashCode(this.array);
	}
	
	public boolean startsWith(ArrayKey<T> key){
		if(key.array.length > array.length){
			return false;
		}
		for(int i = 0; i < key.array.length; i++){
			T e1 = key.array[i];
			T e2 = array[i];
			if(!(e1==null ? e2==null : e1.equals(e2))){
				return false;
			}
		}

		return true;
	}
	
	public String toString(){
		return "["+array[0]+" "+array[1]+"]";
	}
}
