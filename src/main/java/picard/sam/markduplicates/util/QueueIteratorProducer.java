package picard.sam.markduplicates.util;

/*
 * Created by GoodforGod 
 * Date: 24.02.2017
 * Time: 15:02
 */

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.function.Function;

/*
 * DEFAULT COMMENT
 */
public class QueueIteratorProducer<Object, Produced>
{
    private final QueuePeekableIterator<Object> iterator;
    private final Function<QueuePeekableIterator<Object>, Produced> handler;
    private final BlockingQueue<Produced> queue;

    public QueueIteratorProducer(QueuePeekableIterator<Object> iterator, Function<QueuePeekableIterator<Object>, Produced> handler) {
        this.iterator = iterator;
        this.handler = handler;
        this.queue = new LinkedBlockingQueue<>(iterator.CAPACITY);
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
