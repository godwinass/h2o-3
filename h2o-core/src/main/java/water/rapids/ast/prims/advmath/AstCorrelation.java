package water.rapids.ast.prims.advmath;

import water.Key;
import water.MRTask;
import water.fvec.*;
import water.rapids.Env;
import water.rapids.Val;
import water.rapids.vals.ValFrame;
import water.rapids.vals.ValNum;
import water.rapids.ast.AstPrimitive;
import water.rapids.ast.AstRoot;
import water.util.ArrayUtils;
import java.util.Arrays;

/**
 * Calculate Pearson's Correlation Coefficient between columns of a frame
 * <p/>
 * Formula:
 * Pearson's Correlation Coefficient = Cov(X,Y)/sigma(X) * sigma(Y)
 */
public class AstCorrelation extends AstPrimitive {
  @Override
  public String[] args() {
    return new String[]{"ary", "x", "y", "use"};
  }

  private enum Mode {Everything, AllObs, CompleteObs}

  @Override
  public int nargs() {
    return 1 + 3; /* (cor X Y use) */
  }

  @Override
  public String str() {
    return "cor";
  }

  @Override
  public Val apply(Env env, Env.StackHelp stk, AstRoot asts[]) {
    Frame frx = stk.track(asts[1].exec(env)).getFrame();
    Frame fry = stk.track(asts[2].exec(env)).getFrame();
    if (frx.numRows() != fry.numRows())
      throw new IllegalArgumentException("Frames must have the same number of rows, found " + frx.numRows() + " and " + fry.numRows());

    String use = stk.track(asts[3].exec(env)).getStr();
    Mode mode;
    switch (use) {
      case "everything":
        mode = Mode.Everything;
        break;
      case "all.obs":
        mode = Mode.AllObs;
        break;
      case "complete.obs":
        mode = Mode.CompleteObs;
        break;
      default:
        throw new IllegalArgumentException("unknown use mode: " + use);
    }
    return fry.numRows() == 1 ? scalar(frx, fry, mode) : array(frx, fry, mode);
  }

  // Pearson Correlation for one row, which will return a scalar value.
  private ValNum scalar(Frame frx, Frame fry, Mode mode) {
    if (frx.numCols() != fry.numCols())
      throw new IllegalArgumentException("Single rows must have the same number of columns, found " + frx.numCols() + " and " + fry.numCols());
    Vec vecxs[] = frx.vecs();
    Vec vecys[] = fry.vecs();
    double xmean = 0, ymean = 0, xvar = 0, yvar = 0, xsd = 0, ysd = 0, ncols = frx.numCols(), NACount = 0, xval, yval, ss = 0;
    for (int r = 0; r < ncols; r++) {
      xval = vecxs[r].at(0);
      yval = vecys[r].at(0);
      if (Double.isNaN(xval) || Double.isNaN(yval))
        NACount++;
      else {
        xmean += xval;
        ymean += yval;
      }
    }
    xmean /= (ncols - NACount);
    ymean /= (ncols - NACount);

    for (int r = 0; r < ncols; r++) {
      xval = vecxs[r].at(0);
      yval = vecys[r].at(0);
      if (!(Double.isNaN(xval) || Double.isNaN(yval))) {
        //Compute variance of x and y vars
        xvar += Math.pow((vecxs[r].at(0) - xmean), 2);
        yvar += Math.pow((vecys[r].at(0) - ymean), 2);
        //Compute sum of squares of x and y
        ss += (vecxs[r].at(0) - xmean) * (vecys[r].at(0) - ymean);
      }
    }
    xsd = Math.sqrt(xvar / (frx.numRows())); //Sample Standard Deviation
    ysd = Math.sqrt(yvar / (fry.numRows())); //Sample Standard Deviation
    double cor_denom = xsd * ysd;

    if (NACount != 0) {
      if (mode.equals(Mode.AllObs)) throw new IllegalArgumentException("Mode is 'all.obs' but NAs are present");
      if (mode.equals(Mode.Everything)) return new ValNum(Double.NaN);
    }

    for (int r = 0; r < ncols; r++) {
      xval = vecxs[r].at(0);
      yval = vecys[r].at(0);
      if (!(Double.isNaN(xval) || Double.isNaN(yval)))
        ss += (vecxs[r].at(0) - xmean) * (vecys[r].at(0) - ymean);
    }

    return new ValNum(ss / cor_denom); //Pearson's Correlation Coefficient
  }

  // Correlation Matrix.
  // Compute correlation between all columns from each Frame against each other.
  // Return a matrix of correlations which is frx.numCols wide and fry.numCols tall.
  private Val array(Frame frx, Frame fry, Mode mode) {
    Vec[] vecxs = frx.vecs();
    int ncolx = vecxs.length;
    Vec[] vecys = fry.vecs();
    int ncoly = vecys.length;

    if (mode.equals(Mode.AllObs)) {
      for (Vec v : vecxs)
        if (v.naCnt() != 0)
          throw new IllegalArgumentException("Mode is 'all.obs' but NAs are present");
    }

    if(mode.equals(Mode.CompleteObs)){

      CorTaskCompleteObsSum taskCompleteObsMean = new CorTaskCompleteObsSum(ncoly, ncolx).doAll(new Frame(fry).add(frx));
      long NACount = taskCompleteObsMean._NACount;
      double[] ymeans = ArrayUtils.div(taskCompleteObsMean._ysum, fry.numRows() - NACount);
      double[] xmeans = ArrayUtils.div(taskCompleteObsMean._xsum, fry.numRows() - NACount);

      // 1 task with all Xs and Ys
      CorTaskCompleteObs cvs = new CorTaskCompleteObs(ymeans, xmeans).doAll(new Frame(fry).add(frx));

      // 1-col returns scalar 
      if (ncolx == 1 && ncoly == 1) {
        return new ValNum(cvs._cors[0][0] / (cvs._denom[0][0]));
      }

      // Gather all the Xs-vs-Y Correlation arrays; divide by sigma_x and sigma_y
      Vec[] res = new Vec[ncoly];
      Key<Vec>[] keys = Vec.VectorGroup.VG_LEN1.addVecs(ncoly);
      for (int y = 0; y < ncoly; y++)
        res[y] = Vec.makeVec(ArrayUtils.div(cvs.getResult()._cors[y], cvs.getResult()._denom[y]), keys[y]);

      return new ValFrame(new Frame(fry._names, res));

    } else {

      CorTaskEverything[] cvs = new CorTaskEverything[ncoly];

      double[] xmeans = new double[ncolx];
      double[] ymeans = new double[ncoly];
      for (int x = 0; x < ncoly; x++)
        xmeans[x] = vecxs[x].mean();
      for (int y = 0; y < ncoly; y++)
        ymeans[y] = vecys[y].mean();

      // Launch tasks; each does all Xs vs one Y
      for (int y = 0; y < ymeans.length; y++)
        cvs[y] = new CorTaskEverything(ymeans[y], xmeans).dfork(new Frame(vecys[y]).add(frx));

      // 1-col returns scalar
      if (ncolx == 1 && ncoly == 1) {
        return new ValNum(cvs[0].getResult()._cors[0] / cvs[0].getResult()._denom[0]);
      }

      // Gather all the Xs-vs-Y correlation arrays; divide sigma_x and sigma_y
      Vec[] res = new Vec[ncoly];
      Key<Vec>[] keys = Vec.VectorGroup.VG_LEN1.addVecs(ncoly);
      for (int y = 0; y < ncoly; y++)
        res[y] = Vec.makeVec(ArrayUtils.div(cvs[y].getResult()._cors, cvs[y].getResult()._denom), keys[y]);

      return new ValFrame(new Frame(fry._names, res));
    }
  }

  private static class CorTaskEverything extends MRTask<CorTaskEverything> {
    double[] _cors;
    double[] _denom;
    final double _xmeans[], _ymean;

    CorTaskEverything(double ymean, double[] xmeans) {
      _ymean = ymean;
      _xmeans = xmeans;
    }

    @Override
    public void map(Chunk cs[]) {
      final int ncolsx = cs.length - 1;
      final Chunk cy = cs[0];
      final int len = cy._len;
      _cors = new double[ncolsx];
      _denom = new double[ncolsx];
      double sum;
      double varx;
      double vary;
      for (int x = 0; x < ncolsx; x++) {
        sum = 0;
        varx = 0;
        vary = 0;
        final Chunk cx = cs[x + 1];
        final double xmean = _xmeans[x];
        for (int row = 0; row < len; row++) {
          varx += (cx.atd(row) - xmean) * (cx.atd(row) - xmean); //Compute variance for x
          vary += (cy.atd(row) - _ymean) * (cy.atd(row) - _ymean); //Compute variance for y
          sum += (cx.atd(row) - xmean) * (cy.atd(row) - _ymean); //Compute sum of square
        }
        _cors[x] = sum;
        _denom[x] = Math.sqrt(varx) * Math.sqrt(vary);
      }
    }

    @Override
    public void reduce(CorTaskEverything cvt) {
      ArrayUtils.add(_cors, cvt._cors);
      ArrayUtils.add(_denom, cvt._denom);
    }
  }

  private static class CorTaskCompleteObs extends MRTask<CorTaskCompleteObs> {
    double[][] _cors;
    double[][] _denom;
    final double _xmeans[], _ymeans[];

    CorTaskCompleteObs(double[] ymeans, double[] xmeans) {
      _ymeans = ymeans;
      _xmeans = xmeans;
    }

    @Override
    public void map(Chunk cs[]) {
      int ncolx = _xmeans.length;
      int ncoly = _ymeans.length;
      double[] xvals = new double[ncolx];
      double[] yvals = new double[ncoly];
      _cors = new double[ncoly][ncolx];
      _denom = new double[ncoly][ncolx];
      double[] _cors_y;
      double[] sdxy;
      double xval, yval, ymean;
      boolean add;
      int len = cs[0]._len;
      for (int row = 0; row < len; row++) {
        add = true;
        //reset existing arrays to 0 rather than initializing new ones to save on garbage collection
        Arrays.fill(xvals, 0);
        Arrays.fill(yvals, 0);

        for (int y = 0; y < ncoly; y++) {
          final Chunk cy = cs[y];
          yval = cy.atd(row);
          //if any yval along a row is NA, discard the entire row
          if (Double.isNaN(yval)) {
            add = false;
            break;
          }
          yvals[y] = yval;
        }
        if (add) {
          for (int x = 0; x < ncolx; x++) {
            final Chunk cx = cs[x + ncoly];
            xval = cx.atd(row);
            //if any xval along a row is NA, discard the entire row
            if (Double.isNaN(xval)) {
              add = false;
              break;
            }
            xvals[x] = xval;
          }
        }
        //add is true iff row has been traversed and found no NAs among yvals and xvals
        if (add) {
          for (int y = 0; y < ncoly; y++) {
            _cors_y = _cors[y];
            sdxy = _denom[y];
            yval = yvals[y];
            ymean = _ymeans[y];
            for (int x = 0; x < ncolx; x++) {
              _cors_y[x] += (xvals[x] - _xmeans[x]) * (yval - ymean);
              sdxy[x] += Math.sqrt((yval - ymean) * (yval - ymean)) * Math.sqrt((xvals[x] - _xmeans[x]) * (xvals[x] - _xmeans[x]));
            }
          }
        }
      }
    }

    @Override
    public void reduce(CorTaskCompleteObs cvt) {
      ArrayUtils.add(_cors, cvt._cors);
      ArrayUtils.add(_denom, cvt._denom);
    }
  }

  private static class CorTaskCompleteObsSum extends MRTask<CorTaskCompleteObsSum> {
    double[] _xsum, _ysum;
    long _NACount;
    int _ncolx, _ncoly;

    CorTaskCompleteObsSum(int ncoly, int ncolx) {
      _ncolx = ncolx;
      _ncoly = ncoly;
    }

    @Override
    public void map(Chunk cs[]) {
      _xsum = new double[_ncolx];
      _ysum = new double[_ncoly];

      double[] xvals = new double[_ncolx];
      double[] yvals = new double[_ncoly];

      double xval, yval;
      boolean add;
      int len = cs[0]._len;
      for (int row = 0; row < len; row++) {
        add = true;
        //reset existing arrays to 0 rather than initializing new ones to save on garbage collection
        Arrays.fill(xvals, 0);
        Arrays.fill(yvals, 0);

        for (int y = 0; y < _ncoly; y++) {
          final Chunk cy = cs[y];
          yval = cy.atd(row);
          //if any yval along a row is NA, discard the entire row
          if (Double.isNaN(yval)) {
            _NACount++;
            add = false;
            break;
          }
          yvals[y] = yval;
        }
        if (add) {
          for (int x = 0; x < _ncolx; x++) {
            final Chunk cx = cs[x + _ncoly];
            xval = cx.atd(row);
            //if any xval along a row is NA, discard the entire row
            if (Double.isNaN(xval)) {
              _NACount++;
              add = false;
              break;
            }
            xvals[x] = xval;
          }
        }
        //add is true iff row has been traversed and found no NAs among yvals and xvals
        if (add) {
          ArrayUtils.add(_xsum, xvals);
          ArrayUtils.add(_ysum, yvals);
        }
      }
    }

    @Override
    public void reduce(CorTaskCompleteObsSum cvt) {
      ArrayUtils.add(_xsum, cvt._xsum);
      ArrayUtils.add(_ysum, cvt._ysum);
      _NACount += cvt._NACount;
    }
  }

}

