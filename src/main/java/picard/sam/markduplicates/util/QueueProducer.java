package picard.sam.markduplicates.util;

/*
 * Created by GoodforGod 
 * Date: 24.02.2017
 * Time: 15:02
 */

import htsjdk.samtools.util.PeekableIterator;

import java.util.Objects;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.function.Function;

/*
 * DEFAULT COMMENT
 */
public class QueueProducer<Iterable, Produced>
{
    private final int JOB_CAPACITY;

    private final int JOB_CAPACITY_DAFAULT = 8;

    private final PeekableIterator<Iterable> iterator;
    private final Function<PeekableIterator<Iterable>, Produced> handler;
    private final BlockingQueue<Produced> queue;

    public QueueProducer(final PeekableIterator<Iterable> iterator,
                         final Function<PeekableIterator<Iterable>, Produced> handler,
                         int capacity) {
        this.iterator = iterator;
        Objects.requireNonNull(handler, "Handler function must not be null!");
        this.handler = handler;
        this.JOB_CAPACITY = capacity;
        this.queue = new LinkedBlockingQueue<>(JOB_CAPACITY);
        start();
    }

    public QueueProducer(final PeekableIterator<Iterable> iterator,
                         final Function<PeekableIterator<Iterable>, Produced> handler) {
        this.iterator = iterator;
        Objects.requireNonNull(handler, "Handler function must not be null!");
        this.handler = handler;
        this.JOB_CAPACITY = JOB_CAPACITY_DAFAULT;
        this.queue = new LinkedBlockingQueue<>(JOB_CAPACITY);
        start();
    }

    private void start() {
        new Thread(() -> {
            while (iterator.hasNext()) {
                try                             { queue.put(handler.apply(iterator)); }
                catch (InterruptedException e)  { e.printStackTrace(); }
            }
        }).start();
    }

    public boolean hasNext() {
        return queue.peek() != null || iterator.hasNext();
    }

    public Produced peek() {
        return queue.peek();
    }

    public Produced next() {
        try                             { return queue.take(); }
        catch (InterruptedException e)  { e.printStackTrace(); }
        return null;
    }

    public void finish() {
        iterator.close();
    }
}
