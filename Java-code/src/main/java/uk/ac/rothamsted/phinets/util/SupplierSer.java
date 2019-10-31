package uk.ac.rothamsted.phinets.util;

import java.io.Serializable;

import com.google.common.base.Supplier;

public interface SupplierSer<T> extends Supplier<T>, Serializable{

}
