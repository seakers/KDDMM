package seakers.trussaos.initialization;

import java.security.SecureRandom;
import java.util.Random;
import org.apache.commons.math3.random.MersenneTwister;

/**
 * Thread-safe version of Mersenne Twister (used for RandomInitialization)
 * source: https://gist.github.com/dhadka
 */

public class SynchronizedMersenneTwister extends Random {

    private static final long serialVersionUID = -4586969514356530381L;

    private static Random SEEDER;

    private static SynchronizedMersenneTwister INSTANCE;

    private static ThreadLocal<MersenneTwister> LOCAL_RANDOM;

    static {
        SEEDER = new SecureRandom();

        LOCAL_RANDOM = new ThreadLocal<MersenneTwister>() {

            @Override
            protected MersenneTwister initialValue() {
                synchronized (SEEDER) {
                    return new MersenneTwister(SEEDER.nextLong());
                }
            }

        };

        INSTANCE = new SynchronizedMersenneTwister();
    }

    public SynchronizedMersenneTwister() {
        super();
    }

    public static SynchronizedMersenneTwister getInstance() {
        return INSTANCE;
    }

    private MersenneTwister current() {
        return LOCAL_RANDOM.get();
    }

    public synchronized void setSeed(long seed) {
        current().setSeed(seed);
    }

    public void nextBytes(byte[] bytes) {
        current().nextBytes(bytes);
    }

    public int nextInt() {
        return current().nextInt();
    }

    public int nextInt(int n) {
        return current().nextInt(n);
    }

    public long nextLong() {
        return current().nextLong();
    }

    public boolean nextBoolean() {
        return current().nextBoolean();
    }

    public float nextFloat() {
        return current().nextFloat();
    }

    public double nextDouble() {
        return current().nextDouble();
    }

    public synchronized double nextGaussian() {
        return current().nextGaussian();
    }

}