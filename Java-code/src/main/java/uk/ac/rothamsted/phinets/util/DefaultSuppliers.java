package uk.ac.rothamsted.phinets.util;

import java.util.List;
import java.util.Set;

import com.google.common.base.Supplier;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
/**
 * Type-safe Supplier constructors for Guava collections
 * @author lysenkoa
 *
 *		Usage: SetMultimap<String, String> map =  Multimaps.newSetMultimap(Maps.<String, Collection<String>> newHashMap(), DefaultSuppliers.<String>set());
 */
public class DefaultSuppliers {

	public static <T> Supplier<List<T>> list() {
		return new SupplierSer<List<T>>() {
			private static final long serialVersionUID = 1L;

			public List<T> get() {
				return Lists.newArrayList();
			}
		};
	}

	public static <T> Supplier<Set<T>> set() {
		return new SupplierSer<Set<T>>() {
			private static final long serialVersionUID = 1L;

			public Set<T> get() {
				return Sets.newHashSet();
			}
		};
	}
	
	

}
