package picard.sam.markduplicates.util;

/*
 * Created by GoodforGod 
 * Date: 24.02.2017
 * Time: 11:55
 */

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

/*
 * DEFAULT COMMENT
 */
public class QueuePeekableIterator<Object> implements CloseableIterator<Object> {
    private Iterator<Object> iterator;
    private final BlockingQueue<Object> queue;

    private final int QUEUE_CAPACITY_DEFAULT = 8;
    public final int CAPACITY;

    /** Constructs a new iterator that wraps the supplied iterator. */
    public QueuePeekableIterator(Iterator<Object> iterator) {
        this.CAPACITY = QUEUE_CAPACITY_DEFAULT;
        this.queue = new LinkedBlockingQueue<>(CAPACITY);
        this.iterator = iterator;
        //advance();
        start();
    }

    public QueuePeekableIterator(Iterator<Object> iterator, int capacity) {
        this.CAPACITY = capacity;
        this.queue = new LinkedBlockingQueue<>(CAPACITY);
        this.iterator = iterator;
        //advance();
        start();
    }

    /** Closes the underlying iterator. */
    public void close() {
        CloserUtil.close(iterator);
    }

    /** True if there are more items, in which case both next() and peek() will return a value. */
    public boolean hasNext() {
        return this.queue.peek() != null || iterator.hasNext();
    }

    /** Initialize read cycle */
    private void start() {
        new Thread(() -> {
            while (iterator.hasNext()) {
                Object item = advance();

                if(item == null)
                    return;

                try                             { queue.put(item); }
                catch (InterruptedException e)  { e.printStackTrace(); }
            }
        }).start();
    }

    /** Returns the next object and advances the iterator. */
    public Object next() {
        try                             { return queue.take(); }
        catch (InterruptedException e)  { e.printStackTrace(); }
        return null;
    }

    /**
     * Returns the next object but does not advance the iterator. Subsequent calls to peek()
     * and next() will return the same object.
     */
    public Object peek(){
        return queue.peek();
    }

    private Object advance() {
        return this.iterator.hasNext() ? this.iterator.next() : null;
    }

    /** Unsupported Operation. */
    public void remove() {
        throw new UnsupportedOperationException("Not supported: remove");
    }
}
